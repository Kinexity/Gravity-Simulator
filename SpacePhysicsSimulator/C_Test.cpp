#include "C_Test.h"
#include <complex>
#include <valarray>
#include <type_traits>
#include "PureCPPLib/polynomial.h"
#include "PureCPPLib/bit_vector.h"
#include "PureCPPLib/sieve.h"
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/numeric/odeint.hpp>
#include <boost/math/tools/roots.hpp>
#include "matplotlibcpp.h"
#include "gmpxx.h"
namespace pt = matplotlibcpp;

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
	if (mass[1] > mass[0]) {
		std::swap(mass[0], mass[1]);
	}
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
	return std::transform_reduce(std::execution::par, args.begin(), args.end(), static_cast<T>(1.0), std::multiplies<T>(), [&](double _a_) { return (_a_ != arg ? (x - _a_) / (arg - _a_) : static_cast<T>(1.0)); });
}

template<typename T>
T interpolate(const std::vector<double>& args, const std::vector<double>& vals, T x) {
	return std::transform_reduce(std::execution::par, args.begin(), args.end(), vals.begin(), T(), std::plus<T>(), [&](double arg, double val) { return Lagrange_polynomial(args, arg, x) * val;  });
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
	std::generate(std::execution::par, args_interp.begin(), args_interp.end(), interp_arg_gen);
	std::transform(std::execution::par, args_interp.begin(), args_interp.end(), vals_interp.begin(), fn);
	double
		hi = std::max(1.25 * *std::max_element(std::execution::par, vals_interp.begin(), vals_interp.end()), 1.),
		lo = std::min(1.25 * *std::min_element(std::execution::par, vals_interp.begin(), vals_interp.end()), -1.);
	interpolated_polynomial = interpolate(args_interp, vals_interp, polynomial<tp>({ (tp)0.,(tp)1. }));
	std::cout << interpolated_polynomial << '\n';
	std::generate(std::execution::par, draw_args.begin(), draw_args.end(), draw_arg_gen);
	std::transform(std::execution::par, draw_args.begin(), draw_args.end(), vals_draw.begin(), fn);
	std::transform(std::execution::par, draw_args.begin(), draw_args.end(), vals_interp_draw.begin(), [&](auto x) { return std::clamp(interpolate(args_interp, vals_interp, x), lo, hi); });
	std::transform(std::execution::par, draw_args.begin(), draw_args.end(), vals_interp_poly_draw.begin(), [&](auto x) { return (double)std::clamp<tp>(interpolated_polynomial(tp(x)), lo, hi); });
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
	return_value_storing_wrapper fn = ref;
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

void C_Test::run() {
	const size_t N = 1e5;
	//Lagrange_points();
	//interpolation();
	polynomial_speed_test();
	//CVRR_test();
	//zad3();
	//solving_test();
}