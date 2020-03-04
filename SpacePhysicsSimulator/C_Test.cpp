#include "C_Test.h"
#include <complex>
#include <valarray>
#include "matplotlibcpp.h"
using namespace std;
namespace pt = matplotlibcpp;


template <size_t element_index, typename... Types>
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
	return tuple_write<0>(output_stream, data_tuple);
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

class MultiKey3 {
public:
	uint_fast64_t  key[3];

	MultiKey3(uint_fast64_t k1 = -1, uint_fast64_t k2 = -1, uint_fast64_t k3 = -1)
		: key{ k1,k2,k3 } {}

	bool operator<(const MultiKey3& right) const
	{
		if (key[0] == right.key[0]) {
			if (key[1] == right.key[1]) {
				return key[2] < right.key[2];
			}
			else {
				return key[1] < right.key[1];
			}
		}
		else {
			return key[0] < right.key[0];
		}
	}
	bool operator()(const MultiKey3& right) const
	{
		if (key[0] == right.key[0] &&
			key[1] == right.key[1] &&
			key[2] == right.key[2]) return true;
		return false;
	}
};

struct multikey3_less
{
	bool operator()(const MultiKey3& mkey1, const MultiKey3& mkey2) const
	{
		return mkey1.operator<(mkey2); // mkey1 < mkey2;
	}
};

void test(uint_fast64_t thr_ord) {
	static uint_fast64_t i[4];
	static PCL::C_Barrier barrier(4);
	for (i[thr_ord] = 0; i[thr_ord] < 1000000; i[thr_ord]++) {
		barrier.arrive_and_wait();
	}
	std::cout << "watek: " << thr_ord << endl;
	for (uint_fast64_t n = 0; n < 4; n++) {
		std::cout << i[n] << endl;
	}
}

void randomness_test() {
	std::cout << line << std::endl << "randomness_test" << std::endl << line << std::endl;
	std::cout << fixed;
	random_device
		r_d;
	mt19937_64
		gen(r_d());
	uint_fast64_t
		* arr = nullptr,
		counter;
	const uint_fast64_t
		num_of_tests = 10;
	string
		temp;
	for (uint_fast64_t num_of_numbers = 2; num_of_numbers < 50; num_of_numbers++) {
		arr = new uint_fast64_t[num_of_numbers];
		counter = 0;
		for (uint_fast64_t test_ord = 0; test_ord < num_of_tests; test_ord++) {
			for (uint_fast64_t numbers_ord = 0; numbers_ord < num_of_numbers; numbers_ord++) {
				arr[numbers_ord] = gen();
			}
			if (is_sorted(&arr[0], &arr[num_of_numbers - 1])) {
				continue;
			}
			do {
				shuffle(&arr[0], &arr[num_of_numbers - 1], gen);
				counter++;
			} while (!is_sorted(&arr[0], &arr[num_of_numbers - 1]));
		}
		temp = to_string(double(counter) / num_of_tests);
		std::cout << double(counter) / num_of_tests << endl;
	}
	std::cout << scientific;
}

inline double hypotenus(double d1, double d2, double d3) noexcept {
	return sqrt(d1 * d1 + d2 * d2 + d3 * d3);
};

void test_n() {
	thread([&]()->void {
		for (int x = 1; x <= 20; x++) {
			this_thread::sleep_for(1s);
			std::cout << "sekunda " << x << endl;
		}
	}).detach();
	_getch();
	std::cout << "kliknieto klawisz" << endl;
}

void sets_test() {
	PCL::math::interval<double, true, false>
		inter{ 5,67 };
	std::cout << boolalpha << "5 : " << (inter == 5) << endl;
	std::cout << boolalpha << "32 : " << (inter == 32) << endl;
	std::cout << boolalpha << "67 : " << (inter == 67) << endl;
	std::cout << line << endl;
}

void intrinsics_test() {
	std::cout << line << std::endl << "intrinsics_test" << std::endl << line << std::endl;
	uint_fast64_t
		size = 1000000000;
	int_fast64_t
		sum = 0;
	__m256i*
		ptr_avx2 = nullptr;
	__m512i*
		ptr_avx512 = nullptr;
	unique_ptr<int_fast64_t[]>
		array = make_unique<int_fast64_t[]>(size);
	__m256i
		sum_avx2 = { 0,0,0,0 };
	atomic<__m256i>
		atomic_sum_avx2;
	__m512i
		sum_avx512 = _mm512_set1_epi64(int_fast64_t());
	for (uint_fast64_t index = 0; index < size; index++) {
		array[index] = rnd();
	}
	tc().start();
	sum = accumulate(&array[0], &array[size], int_fast64_t());
	tc().stop();
	std::cout << line << endl;
	std::cout << "suma = " << sum << endl;
	std::cout << "accumulate: " << size / tc().measured_timespan().count() << " elem./s" << endl;
	tc().start();
	sum = reduce(execution::par, &array[0], &array[size], int_fast64_t());
	tc().stop();
	std::cout << line << endl;
	std::cout << "suma = " << sum << endl;
	std::cout << "reduce: " << size / tc().measured_timespan().count() << " elem./s" << endl;
	ptr_avx2 = reinterpret_cast<__m256i*>(array.get());
	tc().start();
	sum_avx2 = __m256i{0, 0, 0, 0};
	for (uint_fast64_t mult_index = 0; mult_index < (size / 4); mult_index++) {
		sum_avx2 = _mm256_add_epi64(sum_avx2, ptr_avx2[mult_index]);
	}
	sum = (reinterpret_cast<int_fast64_t*>(&sum_avx2))[0] + (reinterpret_cast<int_fast64_t*>(&sum_avx2))[1] + (reinterpret_cast<int_fast64_t*>(&sum_avx2))[2] + (reinterpret_cast<int_fast64_t*>(&sum_avx2))[3];
	tc().stop();
	std::cout << line << endl;
	std::cout << "suma = " << sum << endl;
	std::cout << "_mm256_add_epi64: " << size / tc().measured_timespan().count() << " elem./s" << endl;
	tc().start();
	sum_avx2 = reduce(execution::par, (__m256i*) & array[0], (__m256i*) & array[size], _mm256_set1_epi64x(0), [&](__m256i v1, __m256i v2) { return _mm256_add_epi64(v1, v2); });
	sum = (reinterpret_cast<int_fast64_t*>(&sum_avx2))[0] + (reinterpret_cast<int_fast64_t*>(&sum_avx2))[1] + (reinterpret_cast<int_fast64_t*>(&sum_avx2))[2] + (reinterpret_cast<int_fast64_t*>(&sum_avx2))[3];
	tc().stop();
	std::cout << line << endl;
	std::cout << "suma = " << sum << endl;
	std::cout << "reduce: " << size / tc().measured_timespan().count() << " elem./s" << endl;
}

void hypot_test() {
	std::cout << line << std::endl << "hypot_test" << std::endl << line << std::endl;
	const uint_fast64_t
		size = 100000000;
	std::vector<double>
		vec;
	double number = double();
	for (uint_fast64_t x = 0; x < size + 9; x++) {
		vec.push_back(static_cast<double>(rnd()));
	}
	tc().start();
	for (uint_fast64_t x = 0; x < size - 5; x++) {
		number = hypot(vec[x] - vec[x + 1], vec[x + 2] - vec[x + 3], vec[x + 4] - vec[x + 5]);
	}
	tc().stop();
	std::cout << "Przepustowosc (hypot): " << size / tc().measured_timespan().count() << " elem./s" << endl;
	tc().start();
	for (uint_fast64_t x = 0; x < size - 5; x++) {
		number = hypot(vec[x] - vec[x + 1], hypot(vec[x + 2] - vec[x + 3], vec[x + 4] - vec[x + 5]));
	}
	tc().stop();
	std::cout << "Przepustowosc (hypot(hypot)): " << size / tc().measured_timespan().count() << " elem./s" << endl;
	tc().start();
	for (uint_fast64_t x = 0; x < size - 5; x++) {
		number = sqrt((vec[x] - vec[x + 1]) * (vec[x] - vec[x + 1]) + (vec[x + 2] - vec[x + 3]) * (vec[x + 2] - vec[x + 3]) + (vec[x + 4] - vec[x + 5]) * (vec[x + 4] - vec[x + 5]));
	}
	tc().stop();
	std::cout << "Przepustowosc (sqrt): " << size / tc().measured_timespan().count() << " elem./s" << endl;
	tc().start();
	for (uint_fast64_t x = 0; x < size - 5; x++) {
		number = sqrt(fma((vec[x] - vec[x + 1]), (vec[x] - vec[x + 1]), fma((vec[x + 2] - vec[x + 3]), (vec[x + 2] - vec[x + 3]), (vec[x + 4] - vec[x + 5]) * (vec[x + 4] - vec[x + 5]))));
	}
	tc().stop();
	std::cout << "Przepustowosc (sqrt(fma)): " << size / tc().measured_timespan().count() << " elem./s" << endl;
	tc().start();
	__m256d
		set_1,
		set_2,
		set_op;
	for (uint_fast64_t x = 0; x < size - 5; x++) {
		set_1 = { vec[x],vec[x + 2],vec[x + 4],double() };
		set_2 = { vec[x + 1],vec[x + 3],vec[x + 5],double() };
		set_op = _mm256_sub_pd(set_2, set_1);
		set_op = _mm256_mul_pd(set_op, set_op);
		number = sqrt(_mm_cvtsd_f64(_mm_add_sd(_mm_add_pd(_mm256_castpd256_pd128(set_op), _mm256_extractf128_pd(set_op, 1)), _mm_unpackhi_pd(_mm_add_pd(_mm256_castpd256_pd128(set_op), _mm256_extractf128_pd(set_op, 1)), _mm_add_pd(_mm256_castpd256_pd128(set_op), _mm256_extractf128_pd(set_op, 1))))));
	}
	tc().stop();
	std::cout << "Przepustowosc (sqrt(intrin)): " << size / tc().measured_timespan().count() << " elem./s" << endl;
	tc().start();
	__m256d
		set[6],
		hypot_set;
	for (uint_fast64_t x = 0; x < size - 5; x++) {
		for (auto index = 0; index < 6; index += 2) {
			set[index] = _mm256_load_pd(&vec[x + index]);
			set[index + 1] = _mm256_load_pd(&vec[x + index + 1]);
			set[index] = _mm256_sub_pd(set[index], set[index + 1]);
		}
		hypot_set = _mm256_hypot_pd(_mm256_hypot_pd(set[0], set[2]), set[4]);
	}
	tc().stop();
	std::cout << "Przepustowosc (hypot(intrin)): " << size / tc().measured_timespan().count() << " elem./s" << endl;
}

void bitset_test() {
	std::cout << line << std::endl << "bitset_test" << std::endl << line << std::endl;
	constexpr size_t
		bitset_size = 200000000;
	size_t
		vb_byte_size = 0;
	vector<bool>
		vector_bitset;
	PCL::bitset<bitset_size>
		my_bitset;
	PCL::C_Time_Counter
		time_counter;
	unique_ptr<bool[]>
		bool_array = make_unique<bool[]>(bitset_size);
	unique_ptr<uint_fast64_t[]>
		index_array = make_unique<uint_fast64_t[]>(bitset_size);
	std::cout << line << endl;
	time_counter.start();
	for (uint_fast64_t index = 0; index < bitset_size; index++) {
		index_array[index] = rnd() % bitset_size;
	}
	time_counter.stop();
	std::cout << "Generacja liczb losowych: " << bitset_size / time_counter.measured_timespan().count() << " elem./s" << endl;
	vector_bitset.reserve(bitset_size);
	std::cout << line << endl;
	std::cout << "PCL::bitset - rozmiar: " << my_bitset.byte_size() << endl;
	time_counter.start();
	for (uint_fast64_t index = 0; index < bitset_size; index++) {
		my_bitset[index] = (index % 2 == 0);
	}
	time_counter.stop();
	std::cout << "PCL::bitset - zapis sekwencyjny: " << bitset_size / time_counter.measured_timespan().count() << " elem./s" << endl;
	time_counter.start();
	for (uint_fast64_t index = 0; index < bitset_size; index++) {
		my_bitset[index_array[index]] = (index % 2 == 0);
	}
	time_counter.stop();
	std::cout << "PCL::bitset - zapis losowy: " << bitset_size / time_counter.measured_timespan().count() << " elem./s" << endl;
	time_counter.start();
	for (uint_fast64_t index = 0; index < bitset_size; index++) {
		bool_array[index] = my_bitset[index];
	}
	time_counter.stop();
	std::cout << "PCL::bitset - odczyt sekwencyjny: " << bitset_size / time_counter.measured_timespan().count() << " elem./s" << endl;
	time_counter.start();
	for (uint_fast64_t index = 0; index < bitset_size; index++) {
		bool_array[index] = my_bitset[index_array[index]];
	}
	time_counter.stop();
	std::cout << "PCL::bitset - odczyt losowy: " << bitset_size / time_counter.measured_timespan().count() << " elem./s" << endl;
	std::cout << line << endl;
	auto ptr_diff = vector_bitset[bitset_size]._Getptr() - vector_bitset[0]._Getptr();
	std::cout << "vector<bool> - rozmiar: " << ptr_diff << endl;
	time_counter.start();
	for (uint_fast64_t index = 0; index < bitset_size; index++) {
		vector_bitset[index] = (index % 2 == 0);
	}
	time_counter.stop();
	std::cout << "vector<bool> - zapis sekwencyjny: " << bitset_size / time_counter.measured_timespan().count() << " elem./s" << endl;
	time_counter.start();
	for (uint_fast64_t index = 0; index < bitset_size; index++) {
		vector_bitset[index_array[index]] = (index % 2 == 0);
	}
	time_counter.stop();
	std::cout << "vector<bool> - zapis losowy: " << bitset_size / time_counter.measured_timespan().count() << " elem./s" << endl;
	time_counter.start();
	for (uint_fast64_t index = 0; index < bitset_size; index++) {
		bool_array[index] = vector_bitset[index];
	}
	time_counter.stop();
	std::cout << "vector<bool> - odczyt sekwencyjny: " << bitset_size / time_counter.measured_timespan().count() << " elem./s" << endl;
	time_counter.start();
	for (uint_fast64_t index = 0; index < bitset_size; index++) {
		bool_array[index] = vector_bitset[index_array[index]];
	}
	time_counter.stop();
	std::cout << "vector<bool> - odczyt losowy: " << bitset_size / time_counter.measured_timespan().count() << " elem./s" << endl;
}

void test() {
	std::cout << line << endl;
	const uint_fast64_t
		size = 10;
	random_device
		r_d;
	mt19937_64
		gen(r_d());
	uint_fast64_t
		a[size];
	double
		n_i_r = 0;
	std::cout << "a[" << size << "] = { ";
	for (uint_fast64_t i = 0; i < size; i++) {
		a[i] = gen() % 10;
		std::cout << a[i] << (i != 9 ? "," : "");
	}
	std::cout << " }" << endl;
	for (uint_fast64_t i = 0; i < size - 1; i++) {
		double l;
		l = (static_cast<double>(a[i]) + a[i + 1]) / 2;
		n_i_r += l;
		std::cout << "(" << a[i] << "+" << a[i + 1] << ")" << "/2" << " = " << l << endl;
		std::cout << "n_i_r = " << n_i_r << endl;
	}
	std::cout << "numerical_integration_r = " << round(n_i_r) << endl;
	//std::cout << "numerical_integration = " << numerical_integration(&a[0], &a[size]) << endl;
}

void speedtest() {
	std::cout << line << std::endl << "speedtest" << std::endl << line << std::endl;
	uint_fast64_t
		arr_size = rnd() % 247 + 3819,
		jump_size = rnd() % 247 + 13,
		position,
		* new_array = new uint_fast64_t[arr_size];
	for (unsigned i = 0; i < arr_size; i++) {
		new_array[i] = 0;
	}
	constexpr uint_fast64_t
		num_of_cycles = 1000000000;
	unique_ptr<uint_fast64_t[]>
		unique_ptr_array = make_unique<uint_fast64_t[]>(arr_size);
	tc().start();
	position = 0;
	for (uint_fast64_t i = 0; i < num_of_cycles; i++) {
		unique_ptr_array[position] += jump_size;
		position += jump_size;
		position %= arr_size;
	}
	tc().stop();
	std::cout << "czas: " << tc().measured_timespan().count() << endl;
	tc().start();
	position = 0;
	for (uint_fast64_t i = 0; i < num_of_cycles; i++) {
		new_array[position] += jump_size;
		position += jump_size;
		position %= arr_size;
	}
	tc().stop();
	std::cout << "czas: " << tc().measured_timespan().count() << endl;
	delete[] new_array;
}

void xml_test(PCL::C_Event_Log_Base& ev_log) {
	std::cout << line << std::endl << "xml_test" << std::endl << line << std::endl;
	string temp;
	PCL::XML_File
		xml(ev_log, filesystem::current_path() / "settings.xml");
	std::cout << xml << endl;
	xml.set_path("trolololo\\inside", true);
	std::cout << "Podaj nazwe node'a so usuniecia: ";
	cin >> temp;
	std::cout << "Node istnieje: " << boolalpha << xml.exists<PCL::XML_Node>(temp) << endl;
	xml.get<>().erase(temp);
	std::cout << "Node istnieje: " << boolalpha << xml.exists<PCL::XML_Node>(temp) << endl;
	xml.set_path("", true);
	std::cout << "Node istnieje: " << boolalpha << xml.exists<PCL::XML_Node>(temp) << endl;
	xml.get<>().erase(temp);
	std::cout << "Node istnieje: " << boolalpha << xml.exists<PCL::XML_Node>(temp) << endl;
	std::cout << xml;
}

void max_test() {
	std::cout << line << std::endl << "max_test" << std::endl << line << std::endl;
	using T = int_fast64_t;
	const size_t
		array_size = 200000000;
	PCL::C_Time_Counter
		time_counter;
	T
		sum = 0;
	unique_ptr<T[]>
		num_array[2] = { make_unique<T[]>(array_size) ,  make_unique<T[]>(array_size) };
	for (uint_fast64_t index = 0; index < array_size; index++) {
		num_array[0][index] = rnd();
		num_array[1][index] = rnd();
	}
	std::cout << line << endl;
	sum = 0;
	time_counter.start();
	for (uint_fast64_t index = 0; index < array_size; index++) {
		sum += std::max(num_array[0][index], num_array[1][index]);
	}
	time_counter.stop();
	std::cout << "Funkcja std::max: " << array_size / time_counter.measured_timespan().count() << " elem./s" << endl;
	std::cout << sum << endl;
	sum = 0;
	time_counter.start();
	for (uint_fast64_t index = 0; index < array_size; index++) {
		sum += (num_array[0][index] > num_array[1][index] ? num_array[0][index] : num_array[1][index]);
	}
	time_counter.stop();
	std::cout << "Operator ternarny: " << array_size / time_counter.measured_timespan().count() << " elem./s" << endl;
	std::cout << sum << endl;
	sum = 0;
	time_counter.start();
	for (uint_fast64_t index = 0; index < array_size; index++) {
		sum += (num_array[0][index] + num_array[1][index] + abs(num_array[0][index] - num_array[1][index]));
	}
	sum /= 2;
	time_counter.stop();
	std::cout << "Clusterfuck: " << array_size / time_counter.measured_timespan().count() << " elem./s" << endl;
	std::cout << sum << endl;
	sum = 0;
	time_counter.start();
	for (uint_fast64_t index = 0; index < array_size; index++) {
		if (num_array[0][index] > num_array[1][index]) {
			sum += num_array[0][index];
		}
		else {
			sum += num_array[1][index];
		}
	}
	sum /= 2;
	time_counter.stop();
	std::cout << "If: " << array_size / time_counter.measured_timespan().count() << " elem./s" << endl;
	std::cout << sum << endl;
}

C_Test::C_Test(PCL::C_Event_Log_Base& event_log_ref) :
	event_log(event_log_ref) {}

inline double pow_sp(double number, uint_fast64_t exponent) noexcept {
	uint_fast64_t max_power_of_2 = floor(log2(exponent));
	double
		result = 1;
	for (uint_fast64_t x = 0; x <= max_power_of_2; x++) {
		if (exponent & (0x1ULL << x)) {
			result *= number;
		}
		number *= number;
	}
	return result;
}

void exponentiation_test() {
	std::cout << line << std::endl << "exp_test" << std::endl << line << std::endl;
	const size_t
		array_size = 2000000;
	unique_ptr<double[]>
		num_array[2] = { make_unique<double[]>(array_size) ,  make_unique<double[]>(array_size) };
	for (auto index = 0; index < array_size; index++) {
		num_array[0][index] = rnd();
	}
	for (auto exponent = 2; exponent <= 20; exponent++) {
		tc().start();
		for (auto index = 0; index < array_size; index++) {
			num_array[1][index] = pow(num_array[0][index], exponent);
		}
		tc().stop();
		std::cout << "pow (" << exponent << "): " << array_size / tc().measured_timespan().count() << " ops/s" << endl;
	}
	for (auto exponent = 2; exponent <= 20; exponent++) {
		tc().start();
		for (auto index = 0; index < array_size; index++) {
			num_array[1][index] = pow_sp(num_array[0][index], exponent);
		}
		tc().stop();
		std::cout << "pow_sp (" << exponent << "): " << array_size / tc().measured_timespan().count() << " ops/s" << endl;
	}
}

union vec_access {
	__m256d avx_vec;
	double array[4];
};

void sum_test() {
	std::cout << line << std::endl << "sum_test" << std::endl << line << std::endl;
	const uint_fast64_t
		size = 200000000;
	double
		sum = 0;
	vec_access
		test;
	__m128d
		temp_vec;
	PCL::C_Time_Counter
		time_counter;
	double
		arr[4] = { double() };
	unique_ptr<double[]>
		elem_array = make_unique<double[]>(size);
	unique_ptr<__m256d[]>
		avx_elem_array = make_unique<__m256d[]>(size / 4);
	for (uint_fast64_t index = 0; index < size; index++) {
		elem_array[index] = rnd();
		if ((index + 1) % 4 == 0) {
			avx_elem_array[(index - 3) / 4] = _mm256_load_pd(&elem_array[index - 3]);
		}
	}
	time_counter.start();
	sum = 0;
	for (uint_fast64_t index = 0; index < size / 4; index++) {
		sum += (reinterpret_cast<double*>(&avx_elem_array[index]))[0] + (reinterpret_cast<double*>(&avx_elem_array[index]))[1] + (reinterpret_cast<double*>(&avx_elem_array[index]))[2] + (reinterpret_cast<double*>(&avx_elem_array[index]))[3];
	}
	time_counter.stop();
	std::cout << "suma = " << sum << endl;
	std::cout << "Przepustowosc (reinterpret_cast): " << size / time_counter.measured_timespan().count() << "elem./s" << endl;
	time_counter.start();
	sum = 0;
	for (uint_fast64_t index = 0; index < size / 4; index++) {
		test.avx_vec = avx_elem_array[index];
		for (auto x = 0; x < 4; x++) {
			sum += test.array[x];
		}
	}
	time_counter.stop();
	std::cout << "suma = " << sum << endl;
	std::cout << "Przepustowosc (vec_access): " << size / time_counter.measured_timespan().count() << "elem./s" << endl;
	time_counter.start();
	sum = 0;
	for (uint_fast64_t index = 0; index < size / 4; index++) {
		temp_vec = _mm256_castpd256_pd128(_mm256_hadd_pd(avx_elem_array[index], avx_elem_array[index]));
		_mm_store_pd(arr, _mm_hadd_pd(temp_vec, temp_vec));
		sum += *arr;
	}
	time_counter.stop();
	std::cout << "suma = " << sum << endl;
	std::cout << "Przepustowosc (_mm256_store_pd): " << size / time_counter.measured_timespan().count() << "elem./s" << endl;
}

void error_test() {
	std::cout << line << std::endl << "error_test" << std::endl << line << std::endl;
	double
		divisor = 8,
		error_magnitude_sum = 0,
		number,
		double_result = 0,
		res_diff = 0;
	double divisor_reciprocal = double(1) / divisor;
	__m256d
		avx_result;
	uint_fast64_t
		error_count = 0;
	const uint_fast64_t
		size = 20000000;
	for (uint_fast64_t counter = 0; counter < size; counter++) {
		number = rnd();
		avx_result = _mm256_mul_pd(
			_mm256_set1_pd(number),
			_mm256_set1_pd(divisor_reciprocal)
		);
		double_result = (number / divisor);
		res_diff = reinterpret_cast<double*>(&avx_result)[0] - double_result;
		if (res_diff != 0) {
			error_count++;
			error_magnitude_sum += log10(abs(res_diff) / double_result);
		}
	}
	std::cout << line << endl;
	std::cout << "Prawidlowosc obliczen: " << static_cast<double>(size - error_count) / size * 100 << " %" << endl;
	std::cout << "Sredni rzad wielkosci bledu obliczen: " << (error_magnitude_sum / size) << endl;
}

void mask_load_test() {
	std::cout << line << std::endl << "mask_load_test" << std::endl << line << std::endl;
	__m128i
		m = _mm_set1_epi64x(0x8000000000000000),
		var = _mm_setr_epi64x(static_cast<long long>(23), static_cast<long long>(45));
	long long
		arr[2] = { 0 ,0 };
	_mm_maskstore_epi64(arr, m, var);
	std::cout << line << endl << arr[0] << "	" << arr[1] << endl;
	std::cout << hex << endl << arr[0] << "	" << arr[1] << dec << endl;
	_mm_store_si128(reinterpret_cast<__m128i*>(arr), var);
	std::cout << line << endl << arr[0] << "	" << arr[1] << endl;
}

double L_radius(const uint_fast64_t L_N, std::array<double, 2> mass, double distance, double accurancy) {
	double
		barycenter_orbit_radius[2] = { double() },
		angular_velocity = double(),
		middle_point = double(),
		edges[2] = { double() };
	constexpr double
		gravitational_constant = 6.6740831313131313131e-11;
	if (mass[1] > mass[0]) {
		swap(mass[0], mass[1]);
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
		default:
			return double();
		}
	};
	auto local_angular_velocity = [&](double local_acceleration, double barycenter_distance) {
		return sqrt(local_acceleration / barycenter_distance);
	};
	switch (L_N) {
	case 1:
		edges[0] = 0;
		edges[1] = distance / (1 + sqrt(mass[1] / mass[0])) - barycenter_orbit_radius[0];
		break;
	case 2:
		edges[0] = barycenter_orbit_radius[1] + double(1e3);
		edges[1] = barycenter_orbit_radius[1] + distance;
		break;
	case 3:
		edges[0] = distance;
		edges[1] = 2 * distance;
		break;
	default:
		return double();
	}
	while (!((edges[1] - edges[0]) < accurancy)) {
		middle_point = (edges[0] + edges[1]) / 2;
		if (local_angular_velocity(local_acceleration(middle_point), middle_point) > angular_velocity) {
			edges[0] = middle_point;
		}
		else {
			edges[1] = middle_point;
		}
	}
	return (edges[0] + edges[1]) / 2 + barycenter_orbit_radius[0];
}

void Lagrange_points() {
	std::array<double, 2>
		mass;
	double
		distance,
		accurancy;
	uint_fast64_t
		L_N;
	std::cout << line << endl;
	std::cout << "Numer punktu: ";
	cin >> L_N;
	std::cout << "Masa 1: ";
	cin >> mass[0];
	std::cout << "Masa 2: ";
	cin >> mass[1];
	std::cout << "Odleglosc: ";
	cin >> distance;
	std::cout << "Dokladnosc: ";
	cin >> accurancy;
	std::cout << "Odleglosc od ciala dominujacego: " << L_radius(L_N, mass, distance, accurancy) << " m" << endl;
}

unsigned char reverse(unsigned char b) {
	b = (b & 0xF0) >> 4 | (b & 0x0F) << 4;
	b = (b & 0xCC) >> 2 | (b & 0x33) << 2;
	b = (b & 0xAA) >> 1 | (b & 0x55) << 1;
	return b;
}

void sum_arrays(size_t exponent) {
	std::cout << line << std::endl << "sum_arrays" << std::endl;
	size_t N = std::pow(10, exponent);
	std::array<std::vector<double>, 3> vecs;
	for (auto& vec : vecs) {
		vec.resize(N);
	}
	for (int i = 0; i < N; i++) {
		vecs[1][i] = rnd();
		vecs[2][i] = rnd();
	}
	tc().start();
	std::transform(std::execution::par, vecs[0].begin(), vecs[0].end(), vecs[2].begin(), [&](double num) {return std::sqrt(num); });
	tc().stop();
	std::cout << "std::transform(std::sqrt): " << N / tc().measured_timespan().count() << "elems/s" << std::endl;
	tc().start();
	for (size_t index = 0; index < N; index += 4) {
		_mm256_store_pd(&vecs[2][index], _mm256_sqrt_pd(_mm256_load_pd(&vecs[0][index])));
	}
	tc().stop();
	std::cout << "_mm256_sqrt_pd: " << N / tc().measured_timespan().count() << "elems/s" << std::endl;
	tc().start();
	std::transform(std::execution::par, (__m256d*) & *vecs[0].begin(), (__m256d*) & *vecs[0].end(), (__m256d*) & *vecs[2].begin(), [&](__m256d v1) { return _mm256_sqrt_pd(v1); });
	tc().stop();
	std::cout << "std::transform(_mm256_sqrt_pd): " << N / tc().measured_timespan().count() << "elems/s" << std::endl;
}

void distance_test() {
	std::cout << line << std::endl << __FUNCTION__ << std::endl << line << std::endl;
	std::set<std::tuple<int, int, int>> s_;
	for (auto i = 0; i < 10; i++) {
		for (auto j = 0; j < 10; j++) {
			for (auto k = 0; k < 10; k++) {
				s_.insert({ i,j,k });
			}
		}
	}
	std::cout << "distance (134,542): " << std::distance(s_.lower_bound({ 1,3,4 }), s_.lower_bound({ 5,4,2 })) << std::endl;
};

void tuple_size_test() {
	std::cout << line << std::endl << __FUNCTION__ << std::endl << line << std::endl;
	std::cout << "sizeof(std::pair<std::tuple<unsigned, int, int>, double>): " << sizeof(std::pair<std::tuple<unsigned, int, int>, double>) << "	" << (sizeof(unsigned) + sizeof(int) + sizeof(int) + sizeof(double)) << std::endl;
	std::cout << "sizeof(std::pair<int, double>): " << sizeof(std::pair<int, double>) << "	" << (sizeof(int) + sizeof(double)) << std::endl;
	std::cout << "sizeof(std::tuple<unsigned, int, int>): " << sizeof(std::tuple<unsigned, int, int>) << "	" << (sizeof(unsigned) + sizeof(int) + sizeof(int)) << std::endl;
}

void matplotlib_test() {
	std::cout << line << std::endl << __FUNCTION__ << std::endl << line << std::endl;
	const int_fast64_t elems = 10000001;
	const double step = 1e-2;
	{
		tc().start();
		std::valarray<double> args, vals;
		args.resize(elems);
		vals.resize(elems);
		std::generate(std::begin(args), std::end(args), [&] {
			static double it = -(elems - 1) / 2 * step;
			return std::exchange(it, it + step);
		});
		vals = std::log(std::sin(std::abs(args)));
		tc().stop();
		std::cout << "valarray: " << elems / tc().measured_timespan().count() << "elems/s" << std::endl;
	}
	tc().start();
	std::vector<double> args, vals;
	args.resize(elems);
	vals.resize(elems);
	std::generate(args.begin(), args.end(), [&] {
		static double it = -(elems - 1) / 2 * step;
		return std::exchange(it, it + step);
	});
	std::transform(args.begin(), args.end(), vals.begin(), [&](double arg) { return std::log(std::sin(std::abs(arg))); });
	tc().stop();
	std::cout << "vector: " << elems / tc().measured_timespan().count() << "elems/s" << std::endl;
	pt::plot(args, vals, "-");
	pt::show();
}

void C_Test::run() {
	/*std::tuple<double, int, std::string> tst_tup;
	std::string str = "342.63 2865 kot_zuje_gume_do_zucia";
	std::stringstream(str) >> tst_tup;
	std::cout << tst_tup << std::endl;
	size_t a;
	cin >> a;*/
	const size_t N = std::pow(10, 7);
	if (false) {
		std::tuple<int, double, unsigned, size_t> t[2] = { {-5, 2.6, 82, 41}, {-5, 2.6, 82, 41} };
		auto res = t[0] + t[1];
		tc().start();
		for (int i = 0; i < N; i++) {
			res += (t[0] + t[1]);
		}
		tc().stop();
		std::cout << std::get<0>(res) << "	" << std::get<1>(res) << "	" << std::get<2>(res) << "	" << std::get<3>(res) << std::endl;
		std::cout << "tuple+: " << N / tc().measured_timespan().count() << " elems/s" << std::endl;
	}
	if (false) {
		uint_fast64_t sum = 0;
		std::map<std::tuple<uint_fast64_t, uint_fast64_t, uint_fast64_t>, uint_fast64_t> mp;
		nested_maps<uint_fast64_t, uint_fast64_t, uint_fast64_t, uint_fast64_t> mp2;
		std::vector<std::array<uint_fast64_t, 4>> arrs;
		arrs.resize(N);
		for (auto& arr : arrs) {
			for (auto& elem : arr) {
				elem = rnd();
			}
		}
		tc().start();
		for (auto& arr : arrs) {
			mp[{arr[0], arr[1], arr[2]}] = arr[3];
		}
		tc().stop();
		std::cout << "array: " << N / tc().measured_timespan().count() << " elems/s" << std::endl;
		tc().start();
		for (auto& arr : arrs) {
			mp2[arr[0]][arr[1]][arr[2]] = arr[3];
		}
		tc().stop();
		std::cout << "nested_maps: " << N / tc().measured_timespan().count() << " elems/s" << std::endl;
		tc().start();
		sum = std::transform_reduce(mp.begin(), mp.end(), uint_fast64_t(), std::plus<>(), [&](auto& elem) {
			return elem.second;
		});
		tc().stop();
		std::cout << "suma - array: " << N / tc().measured_timespan().count() << " elems/s" << std::endl;
		tc().start();
		sum = std::transform_reduce(mp2.begin(), mp2.end(), uint_fast64_t(), std::plus<>(), [](auto& elem)->uint_fast64_t {
			return std::transform_reduce(elem.second.begin(), elem.second.end(), uint_fast64_t(), std::plus<>(), [](auto& elem)->uint_fast64_t {
				return std::transform_reduce(elem.second.begin(), elem.second.end(), uint_fast64_t(), std::plus<>(), [&](auto& elem)->uint_fast64_t {
					return elem.second;
				});
			});
		});
		tc().stop();
		std::cout << "suma - nested_maps: " << N / tc().measured_timespan().count() << " elems/s" << std::endl;
	}
	if (false) {
		uint_fast64_t
			gen_limit = std::numeric_limits<uint_fast64_t>::max(),
			offset = 0;
		std::array<std::vector<uint_fast64_t>, 3> vecs;
		for (auto& v : vecs) {
			v.resize(N);
			for (auto& it : v) {
				it = rnd() % gen_limit + offset;
			}
		}
		std::vector<uint_fast64_t> vec_;
		/*vec_.resize(gen_limit + offset);
		for (auto& elem : vec_) {
			elem = rnd();
		}*/
		std::function ref = std::gcd<uint_fast64_t, uint_fast64_t>;
		C_RV_Storage fn = std::function{ std::gcd<uint_fast64_t, uint_fast64_t> };
		tc().start();
		for (size_t i = 0; i < N; i++) {
			vecs[2][i] = std::gcd(vecs[0][i], vecs[1][i]);
		}
		tc().stop();
		std::cout << "Czas (normalne): " << N / tc().measured_timespan().count() << "elem./s" << endl;
		tc().start();
		for (size_t i = 0; i < N; i++) {
			vecs[2][i] = fn(vecs[0][i], vecs[1][i]);
		}
		tc().stop();
		std::cout << "Czas (CVRR): " << N / tc().measured_timespan().count() << "elem./s" << endl;
	}
	if (false) {
		uint_fast64_t lim;
		bool result = true;
		std::cout << line << endl << "Ile sprawdziæ: ";
		cin >> lim;
		double
			sum = double(1),
			inv_sqrt;
		std::cout << boolalpha;
		for (uint_fast64_t n = 1; n < pow(10, lim); n++) {
			sum *= double(2 * n - 1) / double(2 * n);
			inv_sqrt = 1 / std::sqrt(2 * n);
			std::cout << "L = " << sum << " ; P = " << inv_sqrt << " ; bool(comparison) = " << (result &= (sum < inv_sqrt)) << endl;
		}
		std::cout << line << endl;
		auto fn = [](double x) { return std::sqrt(x - 1) * std::cbrt(log(x)) / std::pow(x, 2. / 3); };
		for (double n__ = double(1); n__ != numeric_limits<double>::infinity(); n__ *= 2) {
			std::cout << scientific << "n = " << n__ << " -> lim = " << fixed << fn(n__) << endl;
		}
		std::cout << std::scientific;
		std::cout << line << std::log(-1. + 0.i) << endl << line << endl;
		//Lagrange_points();
	}
	//xml_test(event_log);
	//bitset_test();
	//tuple_size_test();
	matplotlib_test();
	/*distance_test();
	//sum_arrays(a);
	max_test();
	hypot_test();
	exponentiation_test();
	intrinsics_test();
	sum_test();
	error_test();
	mask_load_test();*/
	std::cout << hex;
	if (false) {
		for (auto x = 0; x < 64; x++) {
			std::cout << "static_cast<bitset_internal>(" << (0x1ULL << x) << ")," << endl;
		}
		std::cout << dec;
		int_fast64_t
			arr[4] = { -23,348290,4198382910, 0 },
			arr_2[4];
		std::cout << line << endl;
		for (auto x = 0; x < 4; x++) {
			std::cout << arr[x] << "	";
		}
		std::cout << endl;
	}
	std::cout << dec;
}