#include "C_Test.h"
#include <complex>
#include <valarray>
#include <type_traits>
#include <numbers>
#include <ranges>
#include "PureCPPLib/polynomial.h"
#include "PureCPPLib/sieve.h"
//#include "particle_container.h"
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/numeric/odeint.hpp>
#include <boost/math/tools/roots.hpp>
#include <Eigen/Dense>
//#include <Eigen/Sparse>
#include <matplot/matplot.h>
#include "PureCPPLib/xoshiro256pp.h"
#include "PureCPPLib/C_thread_set.h"
#undef max
#undef min
//#include <matplotlibcpp.h>
//namespace pt = matplotlibcpp;

static __m256d _mm256_abs_pd(__m256d a) {
	const unsigned long long abs_mask = 0x7FFFFFFFFFFFFFFF;
	const unsigned long long abs_full[8] =
	{ abs_mask, abs_mask, abs_mask, abs_mask, abs_mask, abs_mask, abs_mask,
	   abs_mask };
	return _mm256_and_pd(_mm256_load_pd((double*)abs_full), a);
}

inline
double _mm256_hreduce_pd(__m256d v) {
	__m128d vlow = _mm256_castpd256_pd128(v);
	__m128d vhigh = _mm256_extractf128_pd(v, 1); // high 128
	vlow = _mm_add_pd(vlow, vhigh);     // reduce down to 128
	__m128d high64 = _mm_unpackhi_pd(vlow, vlow);
	return  _mm_cvtsd_f64(_mm_add_sd(vlow, high64));  // reduce to scalar
}

template <typename T> int sgn(T val) {
	return (T() < val) - (val < T());
}

template <size_t element_index = 0, typename... Types>
std::ostream& tuple_write(std::ostream& output_stream, std::tuple<Types...>& data_tuple) {
	if constexpr (element_index < sizeof...(Types)) {
		output_stream << std::get<element_index>(data_tuple) << "	";
		return tuple_write<element_index + 1>(output_stream, data_tuple);
	}
	else {
		return output_stream;
	}
}

template <typename... Types>
std::ostream& operator<<(std::ostream& output_stream, std::tuple<Types...>& data_tuple) {
	return tuple_write<>(output_stream, data_tuple);
}

template <size_t element_index, typename... Types>
std::istream& tuple_read(std::istream& input_stream, std::tuple<Types...>& data_tuple) {
	if constexpr (element_index < sizeof...(Types)) {
		input_stream >> std::get<element_index>(data_tuple);
		return tuple_read<element_index + 1>(input_stream, data_tuple);
	}
	else {
		return input_stream;
	}
}

template <typename... Types>
std::istream& operator>>(std::istream& input_stream, std::tuple<Types...>& data_tuple) {
	return tuple_read<0>(input_stream, data_tuple);
}

C_Test::C_Test(PCL::C_Event_Log_Base& event_log_ref) :
	event_log(event_log_ref) {}

double L_radius(const uint_fast64_t L_N, std::array<double, 2> mass, double distance, double accurancy) {
	double
		barycenter_orbit_radius[2] = { double() },
		angular_velocity = double(),
		middle_point = double(),
		edges[2] = { double() };
	constexpr double
		gravitational_constant = 6.6740831313131313131e-11;
	std::sort(&mass[0], &mass[1]);
	barycenter_orbit_radius[0] = distance * mass[1] / (mass[0] + mass[1]);
	barycenter_orbit_radius[1] = distance * mass[0] / (mass[0] + mass[1]);
	angular_velocity = sqrt(gravitational_constant * (mass[0] + mass[1]) / pow(distance, 3));
	auto local_acceleration = [&](double barycenter_distance) {
		switch (L_N) {
		case 1:
			return gravitational_constant * (mass[0] / pow((barycenter_orbit_radius[0] + barycenter_distance), 2) - mass[1] / pow((barycenter_orbit_radius[1] - barycenter_distance), 2));
		case 2:
			return gravitational_constant * (mass[0] / pow((barycenter_orbit_radius[0] + barycenter_distance), 2) + mass[1] / pow((barycenter_orbit_radius[1] - barycenter_distance), 2));
		case 3:
			return gravitational_constant * (mass[0] / pow((barycenter_orbit_radius[0] - barycenter_distance), 2) + mass[1] / pow((barycenter_orbit_radius[1] + barycenter_distance), 2));
		}
	};
	auto local_angular_velocity = [&](double local_acceleration, double barycenter_distance) {
		return sqrt(local_acceleration / barycenter_distance);
	};
	edges[0] = std::array{ 0., barycenter_orbit_radius[1] + 1e3, distance } [L_N - 1] ;
	edges[1] = std::array{ distance / (1 + sqrt(mass[1] / mass[0])) - barycenter_orbit_radius[0], barycenter_orbit_radius[1] + distance, 2 * distance } [L_N - 1] ;
	auto res = boost::math::tools::bisect([&](double point) { return local_angular_velocity(local_acceleration(middle_point), middle_point) - angular_velocity; },
		edges[0],
		edges[1],
		[&](double min, double max) { return (max - min) < accurancy; });
	return (res.first + res.second) / 2 + barycenter_orbit_radius[0];
}

void Lagrange_points() {
	std::array<double, 2>
		mass;
	double
		distance,
		accurancy;
	uint_fast64_t
		L_N;
	std::cout << line;
	std::cout << "Numer punktu: ";
	std::cin >> L_N;
	if (L_N > 0 && L_N < 4) {
		std::cout << "Masa 1: ";
		std::cin >> mass[0];
		std::cout << "Masa 2: ";
		std::cin >> mass[1];
		std::cout << "Odleglosc: ";
		std::cin >> distance;
		std::cout << "Dokladnosc: ";
		std::cin >> accurancy;
		std::cout << "Odleglosc od ciala dominujacego: " << L_radius(L_N, mass, distance, accurancy) << " m" << '\n';
	}
}

template <typename... Arg_Types>
class differential_equation {
private:
	std::function<double(Arg_Types...)>
		diff_eq_fn;
	const double
		step;
public:
	differential_equation(std::function<double(Arg_Types...)> diff_eq_fn_arg, double step_arg) :
		step(step_arg),
		diff_eq_fn(diff_eq_fn_arg) {}
	~differential_equation() = default;
	double next_value(double current_value, std::tuple<Arg_Types...> fn_args) {
		return current_value + step * std::apply(diff_eq_fn, fn_args);
	}
};

template <size_t N>
using state_type = std::array<double, N>;
template <size_t N>
using error_stepper_type = boost::numeric::odeint::runge_kutta_cash_karp54<state_type<N>>;
template <size_t N>
using controlled_stepper_type = boost::numeric::odeint::controlled_runge_kutta< error_stepper_type<N>>;

class falling_ball_system {
private:
	double
		C = 5.,
		m = 10.,
		g = 10.;
public:
	void operator()(const state_type<2>& x, state_type<2>& dxdt, const double /* t */) {
		dxdt = {
			-g - std::abs(x[0]) * x[0] * C / m,
			x[0]
		};
	}
};

void zad3() {
	tc.start();
	auto height = [&](double t_k) {
		state_type<2> x = { -500. * 1000 / 3600 , 10000. };
		boost::numeric::odeint::integrate_adaptive(controlled_stepper_type<2>(), falling_ball_system(), x, 0., t_k, t_k / 20);
		return x[1];
	};
	auto res = boost::math::tools::bisect(height, 0., 10000., [&](double min, double max) {return max - min < 1e-3; });
	auto time_bisect = std::midpoint(res.first, res.second);
	tc.stop();
	std::cout << "Czas odeint/bisekcja: " << tc.measured_timespan().count() << '\n';
	std::cout << "Flara uderzyla w ziemie po " << time_bisect << " s\n";
}

class var_coeffs_system {
public:
	void operator()(const state_type<3>& x, state_type<3>& dxdt, const double t) {
		dxdt = {
			x[1],
			x[2] * x[0],
			std::cos(t)
		};
	}
};

void solving_test() {
	const auto step = 0.01;
	const auto t_b = -200.;
	const auto t_e = 200.;
	const auto t_0 = 0.;
	std::vector<state_type<3>> states;
	state_type<3> x = { 1., 0., 0. };
	tc.start();
	std::for_each(
		boost::numeric::odeint::make_const_step_iterator_begin(error_stepper_type<3>(), var_coeffs_system(), x, t_0, t_b, -step),
		boost::numeric::odeint::make_const_step_iterator_end(error_stepper_type<3>(), var_coeffs_system(), x),
		[&](auto res) { states.push_back(res); });
	std::reverse(std::execution::par_unseq, states.begin(), states.end());
	states.pop_back();
	x = { 1., 0., 0. };
	std::for_each(
		boost::numeric::odeint::make_const_step_iterator_begin(error_stepper_type<3>(), var_coeffs_system(), x, t_0, t_e, step),
		boost::numeric::odeint::make_const_step_iterator_end(error_stepper_type<3>(), var_coeffs_system(), x),
		[&](auto res) { states.push_back(res); });
	tc.stop();
	std::cout << "Czas odeint: " << tc.measured_timespan().count() << '\n';
	std::fstream file("C:\\Users\\26kuba05\\source\\repos\\SpacePhysicsSimulator\\dane.txt", std::ios::out | std::ios::trunc);
	for (int_fast64_t index = 0; index < states.size(); index++) {
		auto t = t_b + index * (t_e - t_b) / (states.size() - 1);
		if (false) {
			std::cout << "t = " << t << "	y(t) = " << states[index][0] << '\n';
		}
		file << t << "	" << states[index][0] << '\n';
	}
	file.close();
}

template<typename T>
T Lagrange_polynomial(const std::vector<double>& args, double arg, T x) {
	return std::transform_reduce(std::execution::par_unseq, args.begin(), args.end(), static_cast<T>(1.0), std::multiplies<T>(), [&](double _a_) { return (_a_ != arg ? (x - _a_) / (arg - _a_) : static_cast<T>(1.0)); });
}

template<typename T>
T interpolate(const std::vector<double>& args, const std::vector<double>& vals, T x) {
	return std::transform_reduce(std::execution::par_unseq, args.begin(), args.end(), vals.begin(), T(), std::plus<T>(), [&](double arg, double val) { return Lagrange_polynomial(args, arg, x) * val;  });
};

using tp = float256;
double arg_interp, draw_arg;

void interpolation() {
	polynomial<tp>
		interpolated_polynomial;
	const size_t
		N_interp = 101,
		N_draw = 1001;
	const double
		radius = 30.,
		interp_radius = 30.;
	std::vector<double>
		args_interp,
		vals_interp,
		draw_args,
		vals_draw,
		vals_interp_draw,
		vals_interp_poly_draw;
	args_interp.resize(N_interp);
	vals_interp.resize(N_interp);
	draw_args.resize(N_draw);
	vals_draw.resize(N_draw);
	vals_interp_draw.resize(N_draw);
	vals_interp_poly_draw.resize(N_draw);
	arg_interp = -interp_radius;
	draw_arg = -radius;
	auto interp_arg_gen = [&]() {
		static const auto step = 2 * interp_radius / (N_interp - 1);
		return std::exchange(arg_interp, arg_interp + step);
	};
	auto draw_arg_gen = [&]() {
		static const auto step = 2 * radius / (N_draw - 1);
		return std::exchange(draw_arg, draw_arg + step);
	};
	auto fn = [&](auto x) {return std::asinh(x); };
	const double eps = 1e-10;
	auto deriv = [&](auto x) { return (fn(x + eps) - fn(x - eps)) / (2. * eps); };
	std::generate(std::execution::par_unseq, args_interp.begin(), args_interp.end(), interp_arg_gen);
	std::transform(std::execution::par_unseq, args_interp.begin(), args_interp.end(), vals_interp.begin(), fn);
	double
		hi = std::max(1.25 * *std::max_element(std::execution::par_unseq, vals_interp.begin(), vals_interp.end()), 1.),
		lo = std::min(1.25 * *std::min_element(std::execution::par_unseq, vals_interp.begin(), vals_interp.end()), -1.);
	interpolated_polynomial = interpolate(args_interp, vals_interp, polynomial<tp>({ (tp)0.,(tp)1. }));
	std::cout << interpolated_polynomial << '\n';
	std::generate(std::execution::par_unseq, draw_args.begin(), draw_args.end(), draw_arg_gen);
	std::transform(std::execution::par_unseq, draw_args.begin(), draw_args.end(), vals_draw.begin(), fn);
	std::transform(std::execution::par_unseq, draw_args.begin(), draw_args.end(), vals_interp_draw.begin(), [&](auto x) { return std::clamp(interpolate(args_interp, vals_interp, x), lo, hi); });
	std::transform(std::execution::par_unseq, draw_args.begin(), draw_args.end(), vals_interp_poly_draw.begin(), [&](auto x) { return (double)std::clamp<tp>(interpolated_polynomial(tp(x)), lo, hi); });
	//pt::plot(draw_args, vals_draw, "g-");
	//pt::plot(draw_args, vals_interp_draw, "r-");
	//pt::plot(draw_args, vals_interp_poly_draw, "b-");
	//pt::show();
}

void polynomial_speed_test() {
	std::cout << line << __FUNCTION__ << '\n';
	std::vector<polynomial<double>> poly_vec;
	auto real_rnd = std::uniform_real_distribution<double>(1., 2.);
	std::random_device
		r_d;
	std::mt19937_64
		gen{ r_d() };
	const size_t
		N = 1e4,
		max_poly_degree = 10,
		min_poly_degree = 2;
	poly_vec.resize(N);
	std::for_each(std::execution::par_unseq, poly_vec.begin(), poly_vec.end(), [&](polynomial<double>& poly) {
		auto& vec = poly.internal_vector();
		vec.clear();
		auto degree = min_poly_degree + rnd() % (max_poly_degree - min_poly_degree + 1);
		for (auto coef_ind = 0; coef_ind < degree; coef_ind++) {
			vec.push_back(real_rnd(gen));
		}
		});
	auto test_val = real_rnd(gen);
	tc.start();
	auto prod = std::reduce(std::execution::par_unseq, poly_vec.begin(), poly_vec.end(), polynomial<double>({ 1.0 }), std::multiplies<polynomial<double>>());
	tc.stop();
	std::cout << "Czas (standard): " << tc.measured_timespan().count() << '\n';
}

void poly_rem_fn() {
	polynomial<double> poly({ 1.,4.,7.,9. });
	std::cout << "pierwotny: " << poly << '\n';
	std::cout << "pochodna: " << poly.derivate() << '\n';
	std::cout << "calka: " << poly.integral() << '\n';
	polynomial<double> poly_2({ -6.,-5.,2.,1. });
	auto steps = 1000;
	std::fstream file;
	file.open("C:\\Users\\26kuba05\\source\\repos\\SpacePhysicsSimulator\\dane.txt", std::ios::out | std::ios::trunc);
	for (auto i = 0; i <= steps; i++) {
		auto x = -10. + i * 20. / steps;
		polynomial<double> poly_3({ x, 1., 2. });
		auto [res, rem] = poly_2 / poly_3;
		file << x << "	" << rem[0] << '\n';
	}
	file.close();
}

uint_fast64_t binomial_coefficient(uint_fast64_t n, uint_fast64_t k) {
	[[unlikely]] if (k > n) {
		return 0;
	}
	if (n - k < k) {
		k = n - k;
	}
	uint_fast64_t result = 1;
	for (uint_fast64_t i = 1; i <= k; i++) {
		result *= (n - i);
		result /= i;
	}
	if (n != k) {
		result /= (n - k);
		result *= n;
	}
	return result;
};

void CVRR_test() {
	uint_fast64_t
		gen_limit = 1000,
		offset = 0;
	const size_t
		N = 1e7;
	std::array<std::vector<uint_fast64_t>, 2> vecs;
	std::vector<uint_fast64_t> res;
	res.resize(N);
	for (auto& v : vecs) {
		v.resize(N);
	}
	std::for_each(std::execution::par_unseq, res.begin(), res.end(), [&](auto& elem) {
		auto i = std::distance(&*res.begin(), &elem);
		std::tie(vecs[1][i], vecs[0][i]) = std::minmax(rnd() % gen_limit + offset, rnd() % gen_limit + offset);
		});
	std::function ref = binomial_coefficient;
	return_value_wrapper fn = ref;
	tc.start();
	for (size_t i = 0; i < N; i++) {
		res[i] = binomial_coefficient(vecs[0][i], vecs[1][i]);
	}
	tc.stop();
	std::cout << "Czas (normalne): " << N / tc.measured_timespan().count() << "elem./s" << '\n';
	tc.start();
	for (size_t i = 0; i < N; i++) {
		res[i] = fn(vecs[0][i], vecs[1][i]);
	}
	tc.stop();
	std::cout << "Czas (CVRR): " << N / tc.measured_timespan().count() << "elem./s" << '\n';
}

void intrin_test() {
	std::cout << line << __FUNCTION__ << '\n';
	const size_t
		N = 1e8;
	std::array<std::vector<int_fast32_t>, 2> vecs;
	std::vector<double> res;
	res.resize(N);
	for (auto& v : vecs) {
		v.resize(N);
		std::for_each(std::execution::par_unseq, v.begin(), v.end(), [&](auto& elem) {
			elem = (rnd() % 100 + 1);// 2.;
			});
	}
	//tc.start();
	//for (auto i = 0; i < N; i += 4) {
	//	auto tmp = _mm256_load_pd(&vecs[0][i]);
	//	auto tmp2 = _mm256_load_pd(&vecs[1][i]);
	//	_mm256_store_pd(&res[i], _mm256_add_pd(tmp, tmp2));
	//}
	//tc.stop();
	//std::cout << "_mm256_add_pd: " << N / tc.measured_timespan().count() << "elem/s\n";
	//tc.start();
	//for (auto i = 0; i < N; i += 4) {
	//	auto tmp = _mm256_load_pd(&vecs[0][i]);
	//	auto tmp2 = _mm256_load_pd(&vecs[1][i]);
	//	_mm256_store_pd(&res[i], _mm256_mul_pd(tmp, tmp2));
	//}
	//tc.stop();
	//std::cout << "_mm256_mul_pd: " << N / tc.measured_timespan().count() << "elem/s\n";
	//tc.start();
	//for (auto i = 0; i < N; i += 4) {
	//	auto tmp = _mm256_load_pd(&vecs[0][i]);
	//	auto tmp2 = _mm256_load_pd(&vecs[1][i]);
	//	_mm256_store_pd(&res[i], _mm256_div_pd(tmp, tmp2));
	//}
	//tc.stop();
	//std::cout << "_mm256_div_pd: " << N / tc.measured_timespan().count() << "elem/s\n";
	tc.start();
	for (auto i = 0; i < N; i++) {
		vecs[0][i] = vecs[1][i] * 2;
	}
	tc.stop();
	std::cout << "pre: " << N / tc.measured_timespan().count() << "elem/s\n";
	std::cout << "checksum: " << std::reduce(std::execution::par_unseq, vecs[0].begin(), vecs[0].end(), 0, std::plus<>()) << '\n';
	tc.start();
	for (auto i = 0; i < N; i += 8) {
		_mm256_store_si256((__m256i*) & vecs[0][i], _mm256_mul_epi32(_mm256_load_si256((__m256i*) & vecs[1][i]), _mm256_set1_epi32(2)));
	}
	tc.stop();
	std::cout << "post: " << N / tc.measured_timespan().count() << "elem/s\n";
	std::cout << "checksum: " << std::reduce(std::execution::par_unseq, vecs[0].begin(), vecs[0].end(), 0, std::plus<>()) << '\n';
}

class particle_system {
private:
	double
		C = 5.,
		g = 10.;
public:
	particle_system(double C_arg = 5., double g_arg = 10.) : C(C_arg), g(g_arg) {};
	void operator()(const state_type<2>& x, state_type<2>& dxdt, const double /* t */) {
		dxdt = {
			-g - x[0] * x[0] * C,
			x[0]
		};
	}
};

void solution_test() {
	std::vector<double>
		h_result,
		time_points;
	auto height = [&](double t_k) {
		state_type<2> x = { 10 , 10000. };
		boost::numeric::odeint::integrate_adaptive(controlled_stepper_type<2>(), particle_system(), x, 0., t_k, t_k / 20);
		return x[1];
	};
	auto res = boost::math::tools::bisect(height, 0., 10000., [&](double min, double max) {return max - min < 1e-3; });
	auto time_bisect = std::midpoint(res.first, res.second);
	const auto steps = 1001;
	for (auto step = 0; step < steps; step++) {
		auto t = time_bisect * step / (steps - 1);
		time_points.push_back(t);
		h_result.push_back(height(t));
	}
	std::cout << "Czasteczka uderzyla w ziemie po " << time_bisect << " s\n";
	//matplot::plot(time_points, h_result, "r-");
	//matplot::show();
}

void rename_files() {
	std::vector<std::pair<std::wstring, std::wstring>> vec;
	std::filesystem::path dir("D:\\YTBUP\\MxR Mods");
	for (auto& file : std::filesystem::directory_iterator(dir)) {
		auto old_fname = file.path().filename().wstring();
		auto num_str = old_fname.substr(0, 4);
		auto name_str = old_fname.substr(4);
		int num = std::stoi(num_str) - 373;
		auto new_fname = (std::wstringstream() << std::setw(4) << std::setfill(L'0') << num << L'.').str();
		new_fname += name_str;
		vec.push_back({ old_fname, new_fname });
	}
	for (auto& [old_fname, new_fname] : vec) {
		std::filesystem::rename(dir / old_fname, dir / new_fname);
	}
}

void filter_data() {
	std::filesystem::path dir("C:\\Users\\26kuba05\\source\\NewFolder1");
	for (auto& file : std::filesystem::directory_iterator(dir)) {
		std::vector<std::string> lines;
		auto path = file.path();
		std::fstream file;
		file.open(path, std::ios::in);
		std::string line;
		while (file >> line) {
			lines.push_back(line);
		}
		file.close();
		std::vector<std::string> new_lines{ lines.begin() + 2,lines.end() };
		file.open(path, std::ios::out | std::ios::trunc);
		for (auto& elem : new_lines) {
			std::replace(elem.begin(), elem.end(), ',', ' ');
			file << elem << '\n';
		}
		file.close();
	}
}

class ChebyshevInterpolant {
public:
	// do not allow default constructor
	ChebyshevInterpolant() = delete;
	ChebyshevInterpolant(std::function<double(double)> f, std::pair<double, double> section, int degree);
	double Evaluate(double x);
private:
	// section represent begining and end of section on which we interpolate
	std::pair<double, double> section;
	// maximal degree of polynomial
	int degree;
	// pointer to the function being interpolated
	std::function<double(double)> f;
	// list of c_k coefficients
	std::vector<double> c;
	polynomial<double> interpolation_polynomial;
};

ChebyshevInterpolant::ChebyshevInterpolant(std::function<double(double)> f, std::pair<double, double> section, int degree) {
	double xj, xxj, val;
	for (int i = 0; i <= degree; i++) {
		double ci = 0.0;
		for (int j = 0; j < degree; j++) {
			xj = (std::numbers::pi * (j + 0.5)) / degree;
			xxj = 0.5 * (section.first + section.second) + 0.5 * (section.second - section.first) * cos(xj);
			val = f(xxj);
			xxj = cos(i * xj);
			ci = ci + val * xxj;
		}
		ci = (ci * 2.0) / degree;
		c.push_back(ci);
	}
	polynomial<double>
		T_n_1(1, 0),
		T_n(1, 1);
	interpolation_polynomial += c[1] * T_n + c[0] * T_n_1 - 0.5 * c[0];
	for (auto& k : c) {
		auto i = std::distance(&*c.begin(), &k);
		std::cout << "c_" << i << "= " << k << std::endl;
		if (i >= 2) {
			T_n_1 = std::exchange(T_n, 2 * polynomial<double>(1, 1) * T_n - T_n_1);
			interpolation_polynomial += T_n * c[i];
		}
	}
}

double ChebyshevInterpolant::Evaluate(double x) {
	return interpolation_polynomial(x);
};

int czeb_test() {
	double a = 0.5;
	double b = 5.0;
	auto p = std::make_pair(a, b);
	int degree = 4;
	ChebyshevInterpolant Inter1([&](double x) { return std::log(x); }, p, degree);
	double x = 2.5;
	double FunVal = Inter1.Evaluate(x);
	std::cout << "Moja wartosc szacowana = ln(" << x << ") =" << FunVal << std::endl;
	std::cout << "Porownanie z wartoscia z math.h = ln(" << x << ") = " << log(x) << std::endl;
	std::vector<std::pair<double, double>> proba;
	double odcinek_prob = (b - a) / 1000;
	for (int i = 1; i <= 1000; i++) {
		double xi = a + i * odcinek_prob;
		double fxi = Inter1.Evaluate(xi);
		auto p1 = std::make_pair(xi, fxi);
		proba.push_back(p1);
	}
	double max_niep = 0.0;
	for (int e = 1; e < proba.size(); e++) {
		double x1 = proba.at(e).second - log(proba.at(e).first);
		double x2 = proba.at(e - 1).second - log(proba.at(e - 1).first);
		if (abs(x1) > abs(x2)) { max_niep = abs(x1); }
		else { max_niep = abs(x2); }
	}
	std::cout << "maksimum bledu aproksymacji = " << max_niep << std::endl;
	return 0;

}

double fun(int i) {
	return 1 / (double)i / (double)i;
}

std::vector<double> Sum_series(std::function<double(int)> f, int i) {
	std::vector<double> v(i), v1(i);
	std::ranges::for_each(v1, [&](auto& x) { x = f(std::distance(&*v1.begin(), &x) + 1); });
	std::partial_sum(v1.begin(), v1.end(), v.begin());
	return v;
}

int silnia(int k) {
	int sil = 1;
	while (k > 1) {
		sil = sil * k;
		k--;
	}
	return sil;
}

double Rich_extrapolation(std::vector<double> v, int k) {
	double R_k = 0.0;
	int N = v.size() - k;
	for (int i = 0; i < k; i++) {
		R_k += pow((N + i), i) * pow(-1, (k + i)) * v.at(N + i - 1) / (silnia(i) * silnia(k - i));
	}
	return R_k;
}

int rich_test() {
	double porownanie = std::pow(std::numbers::pi, 2) / 6;
	int N = 1000;
	std::vector<double> Sum_vec = Sum_series(fun, N);
	double R1 = N * Sum_vec.at(N - 1) - (N - 1) * Sum_vec.at(N - 2);
	double R2 = (N * N * Sum_vec.at(N - 1) - 2 * (N - 1) * (N - 1) * Sum_vec.at(N - 2) + (N - 2) * (N - 2) * Sum_vec.at(N - 3)) / 2;
	std::cout << "(Na piechote) R1-pi^2/6 = " << R1 - porownanie << std::endl;
	std::cout << "(Na piechote) R2-pi^2/6 = " << R2 - porownanie << std::endl;
	std::cout << "R1-pi^2/6 = " << Rich_extrapolation(Sum_vec, 1) - porownanie << std::endl;
	std::cout << "R2-pi^2/6 = " << Rich_extrapolation(Sum_vec, 2) - porownanie << std::endl;

	return 0;
}

void sieve_test() {
	std::cout << line;
	std::cout << std::fixed;
	unsigned int
		i_start = 1,
		i_stop = 1;
	std::cin >> i_start >> i_stop;
	for (auto i = i_start; i < i_stop; i++) {
		tc.start();
		const auto&& res = PCL::sieve(std::pow(10, i));
		tc.stop();
		std::cout << i << '\t' << tc.measured_timespan().count() << " s \n";
	}
	std::cout << std::scientific;
}

//void mp_test() {
//	std::vector<double> vec(10);
//	std::iota(vec.begin(), vec.end(), 0.);
//	matplot::plot(vec, "g-");
//	matplot::show();
//}

void print_img(const std::vector<double>& img) {
	for (int i = 0; i < 28; i++) {
		for (int j = 0; j < 28; j++) {
			auto px = img[i * 28 + j];
			std::cout << (px > 0.5 ? "#" : (px > 0 ? "*" : " "));
		}
		std::cout << '\n';
	}
}

void get_libs_names() {
	std::cout << std::filesystem::current_path().string() << '\n';
	std::fstream file;
	file.open("libs.txt", std::ios::in);
	std::string line, lib_name;
	while (std::getline(file, line)) {
		std::stringstream(line) >> lib_name;
		std::cout << lib_name << "	";
	}
	std::cout << '\n';
}

void prob_test() {
	std::uniform_real_distribution<double> uni{ -4.,4. };
	size_t N = 1e8;
	std::vector<double> a(N);
	tc.start();
	PCL::C_thread_set ts;
	auto threads = std::thread::hardware_concurrency();
	std::atomic<int> gtz = 0;
	std::mutex mtx;
	std::function fn = [&]() {
		xoshiro256pp gen;
		int local_sum = 0;
		for (int i = 0; i < N / threads; i++) {
			auto b = uni(gen);
			local_sum += b * b - 4 * uni(gen) * uni(gen) < 0.;
		}
		gtz += local_sum;
	};
	ts.run(fn, threads);
	tc.stop();
	std::cout << tc.measured_timespan().count() << '\t' << double(gtz) / N << '\n';
};

void ts_test() {
	auto ftr = [&](size_t n) {
		double f = 1.;
		for (size_t i = 1; i < n; i++) {
			f *= i;
		}
		return f;
	};
	tc.start();
	size_t N = 13;
	const std::uniform_real_distribution<double> uni{ 0.,1. };
	std::vector<double> x, y, path_lenghts;
	//std::vector<std::vector<double>> dist;
	Eigen::MatrixXd dist = Eigen::MatrixXd::Zero(N, N);
	std::vector<size_t> order(N);
	for (size_t i = 0; i < N; i++) {
		x.push_back(uni(rnd));
		y.push_back(uni(rnd));
	}
	for (size_t i = 0; i < N; i++) {
		for (size_t j = 0; j < N; j++) {
			dist(i, j) = std::hypot(x[i] - x[j], y[i] - y[j]);
		}
	}
	std::iota(order.begin(), order.end(), 0);
	double lenght, min_lenght = std::numeric_limits<double>::infinity();
	do {
		lenght = 0.;
		for (size_t i = 0; i < N; i++) {
			lenght += dist(order[i], order[(i + 1) % N]);
		}
		path_lenghts.push_back(lenght);
	} while (std::next_permutation(order.begin() + 1, order.end()));
	tc.stop();
	double avg = std::reduce(std::execution::par_unseq, path_lenghts.begin(), path_lenghts.end(), 0., std::plus<>()) / ftr(N);
	std::sort(std::execution::par_unseq, path_lenghts.begin(), path_lenghts.end());
	std::cout << "Time: " << tc.measured_timespan().count() << " s\n";
	std::cout << "Average: " << avg << "\n";
	std::cout << "Min: " << *std::min_element(std::execution::par_unseq, path_lenghts.begin(), path_lenghts.end()) << "\n";
};

void basis_test() {
	// Define the matrix A and the vector b
	Eigen::MatrixXd A(3, 3);
	A << 1, 2, 3, 4, 5, 6, 7, 8, 9;
	Eigen::VectorXd b(3);
	b << 1, 2, 3;

	// Solve the system of linear equations using the QR decomposition method
	Eigen::VectorXd x = A.colPivHouseholderQr().solve(b);

	// Find the null space of A using the fullPivHouseholderQr() method
	Eigen::FullPivHouseholderQR<Eigen::MatrixXd> qr(A);
	Eigen::MatrixXd Q = qr.matrixQ();
	Eigen::MatrixXd nullspace = Q.rightCols(A.cols() - qr.rank());

	// Print the basis vectors of the null space
	std::cout << "Basis vectors of the null space:\n" << nullspace << std::endl;
}

void solve_test() {
	Eigen::Matrix3f m;
	Eigen::Vector3f v;
	m << 1, 0, 0, 0, 1, 1, 0, 1, 1;
	v << 1, 1, 1;
	std::cout << m << "\n\n" << v << "\n\n";
	auto a = m.llt().solve(v);
	std::cout << a << "\n\n";
};

// Define constants
const double f1 = 87.98;
const double sigma_f1 = 0.16;
const double f2 = 16.36;
const double sigma_f2 = 0.07;

std::vector<bool> indices;

// MC how many muons pass both layers
double MC_count_pass(int N, double a, double h, double gamma) {
	std::uniform_real_distribution<double> dis(0.0, 1.0);
	double prob_sum = 0;
	std::map<std::thread::id, xoshiro256pp> gens;
	auto h_a_ratio = h / a;

	prob_sum = std::transform_reduce(std::execution::par, indices.begin(), indices.end(), double(), std::plus<double>(), [&](bool ph)->double {
		auto& gen = gens[std::this_thread::get_id()];
		double inv_cos_2_theta = std::pow(dis(gen), - 2 / (gamma + 1));
		auto _tan_ = h_a_ratio * std::sqrt(inv_cos_2_theta - 1);
		double phi_ = dis(gen) * std::numbers::pi / 2;
		double x_b = 1 - _tan_ * cos(phi_);
		double y_b = 1 - _tan_ * sin(phi_);

		return (x_b >= 0) * (y_b >= 0) * x_b * y_b;
		});
	return prob_sum / N;
}

// Calculate ratio
double ratio(double gamma, uint64_t N = 10e7) {
	return MC_count_pass(N, 1.0, 0.5, gamma) / MC_count_pass(N, 1.0, 2.0, gamma);
}

void MC_prob3() {
	PCL::C_Time_Counter tc;
	tc.start();

	double f_ratio = f1 / f2;
	double gamma_init = 1.0; // Initial guess for gamma
	double tol = 1e-5; // Tolerance for convergence
	auto N = 10e6;
	indices.resize(N);

	auto result = boost::math::tools::bisect([&](double gamma) { return ratio(gamma, N) - f_ratio; },
		0.1,
		5.,
		[&](double min, double max) { return (max - min) < tol; });
	
	std::cout << "Gamma from MC simulation: " << result.first << std::endl;
	auto phi_tot = (1 / MC_count_pass(N, 1.0, 0.5, result.first)) * f1;
	std::cout << "Phi_tot from MC: " << phi_tot << std::endl;

	tc.stop();

	std::cout << tc.measured_timespan().count() << " s\n";

	for (double gamma = 1.7; gamma < 1.8; gamma += 0.01) {
		std::cout << gamma << "\t\t" << ratio(gamma, N) - f_ratio << "\n";
	}
}


void C_Test::run() {
	//Lagrange_points();
	//interpolation();
	//polynomial_speed_test();
	//solution_test();
	//intrin_test();
	//additional_CPU_info();
	//CVRR_test();
	//zad3();
	//solving_test();
	//filter_data();
	//rich_test();
	//mp_test();
	//prob_test();
	//ts_test();
	//basis_test();
	//solve_test();
	MC_prob3();
}

class Minesweeper_field {
private:
	Eigen::MatrixXi
		bomb_field,
		known_field;
	Eigen::MatrixXf
		bomb_probability_field;
	//Eigen::SparseMatrix<int> sp_mat;
	int32_t
		rows,
		cols,
		mine_count;
	xoshiro256pp
		rng;
	bool
		game_lost;
	bool within_bounds(int32_t row, int32_t col);
	void set_mines(int32_t row, int32_t col);
	void check_mine_equality(int32_t row, int32_t col);
	void explore(int32_t row, int32_t col);
public:
	Minesweeper_field() = default;
	~Minesweeper_field() {};
	void initialize_field();
	bool solve_field();
};

bool Minesweeper_field::within_bounds(int32_t row, int32_t col) {
	return std::clamp(row, 0, rows) == row && std::clamp(col, 0, cols) == col;
}

void Minesweeper_field::set_mines(int32_t row, int32_t col) {
	for (int row_with_shift = std::max(row - 1, 0); row_with_shift < std::min(row + 2, rows); row_with_shift++) {
		for (int col_with_shift = std::max(col - 1, 0); col_with_shift < std::min(col + 2, rows); col_with_shift++) {
			known_field(row_with_shift, col_with_shift) -= (known_field(row_with_shift, col_with_shift) > 0) || (known_field(row_with_shift, col_with_shift) == -1);
		}
	}
}

void Minesweeper_field::check_mine_equality(int32_t row, int32_t col) {
	if (known_field(row, col) > 0) {
		int count = 0; //unknown fields counter
		for (int row_with_shift = std::max(row - 1, 0); row_with_shift < std::min(row + 2, rows); row_with_shift++) {
			for (int col_with_shift = std::max(col - 1, 0); col_with_shift < std::min(col + 2, rows); col_with_shift++) {
				count += (col_with_shift != col || row_with_shift != row) && (known_field(row_with_shift, col_with_shift) == -1);
				mine_count -= (known_field(row_with_shift, col_with_shift) == -1);
			}
		}
		if (count == known_field(row, col)) {
			set_mines(row, col);
		}
		for (int row_with_shift = std::max(row - 1, 0); row_with_shift < std::min(row + 2, rows); row_with_shift++) {
			for (int col_with_shift = std::max(col - 1, 0); col_with_shift < std::min(col + 2, rows); col_with_shift++) {
				check_mine_equality(row_with_shift, col_with_shift);
			}
		}
	}
	else if (known_field(row, col) == 0) {
		explore(row, col);
	}
}

void Minesweeper_field::explore(int32_t row, int32_t col) {

}

void Minesweeper_field::initialize_field() {
	bomb_field = Eigen::MatrixXi::Zero(rows, cols);
	known_field = -Eigen::MatrixXi::Ones(rows, cols);
	std::vector<size_t>
		indices(rows * cols),
		sampled_indices(mine_count);
	std::iota(indices.begin(), indices.end(), 0);
	std::sample(indices.begin(), indices.end(), sampled_indices.begin(), mine_count, rng);
	for (auto i : sampled_indices) {
		bomb_field(i) = 1;
	}
}

bool Minesweeper_field::solve_field() {
	do {
		//iterate over known fields to check if all unknown fields around each could be bombs
		for (int row = 0; row < rows; row++) {
			for (int col = 0; col < cols; col++) {
				if (known_field(row, col) > 0) {
					check_mine_equality(row, col);
				}
			}
		}
	} while (mine_count > 0 && !game_lost);
	return false;
}
