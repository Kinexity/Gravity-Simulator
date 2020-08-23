#pragma once
#ifndef C_obj
#define C_obj
#define _ENABLE_EXTENDED_ALIGNED_STORAGE
#include <thread>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <cstdint>
#include <array>
#include <PureCPPLib/C_Random.h>
#undef max

template <size_t dims>
class C_Object {
private:
	std::array<double, dims>
		position,
		velocity;
	double
		mass = double(),
		radius = double(),
		standard_grav_param = double();
public:
	uint_fast64_t
		object_id = 0;
	C_Object() = default;
	C_Object(C_Object<dims>& obj);
	C_Object(C_Object<dims>&&) = default;
	C_Object<dims>& operator=(const C_Object<dims>&) = default;
	~C_Object() = default;
	bool
		operator!=(C_Object<dims>& obj),
		operator==(C_Object<dims>& obj);
	void
		reset(),
		data_rand(uint_fast64_t noo),
		calculate_sgp();
	double
		operator()(C_Object<dims>& obj_2),
		velocity_val();
	template <size_t dimensions>
	friend class
		C_Universe;
	friend class
		C_Analysis;
	friend std::ostream& operator<<(std::ostream& str, C_Object<dims>& obj);
	friend std::istream& operator>>(std::istream& str, C_Object<dims>& obj);
};

template <size_t dims>
inline void C_Object<dims>::calculate_sgp() {
	standard_grav_param = mass * 6.67430151515151515151515e-11;
}

template <size_t dims>
inline double C_Object<dims>::operator()(C_Object<dims>& obj_2) {
	return std::sqrt(std::transform_reduce(std::execution::seq, position.begin(), position.end(), obj_2.position.begin(), 0.0, std::plus<>(), [&](auto p1, auto p2) {return std::pow(p1 - p2, 2); }));
}

template <size_t dims>
inline double C_Object<dims>::velocity_val() {
	return std::sqrt(std::transform_reduce(std::execution::seq, velocity.begin(), velocity.end(), 0.0, std::plus<>(), [&](auto p1) {return std::pow(p1, 2); }));
}

template <size_t dims>
inline void C_Object<dims>::data_rand(uint_fast64_t noo) {
	static std::random_device rd;
	static std::uniform_real_distribution urdist(0.19884, std::nextafter(1.9884, std::numeric_limits<double>::max()));
	mass = urdist(rd) * pow(10, 29);
	radius = pow(10, 1);
	auto x = static_cast<uint_fast64_t>(sqrt(cbrt(noo)));
	auto min_pow = 6;
	for (auto& elem : position) {
		elem = ((static_cast<double>(rnd() % 9000000 + 1000000)) * pow(10, rnd() % 5 + min_pow) * ((rnd() % 2 == 1) ? 1.0 : -1.0));
	}
	for (auto& elem : velocity) {
		elem = (static_cast<double>(rnd() % 9000 + 1000) * ((rnd() % 2 == 1) ? 1.0 : -1.0));
	}
}

//template <size_t dims>
//C_Object<dims>::C_Object(C_Universe& universe_ref) :
//parent_universe(universe_ref) {}

template <size_t dims>
inline C_Object<dims>::C_Object(C_Object<dims>& obj) {
	memcpy(&position, &obj.position, ((int)&buffer - (int)&position));
}

template <size_t dims>
inline bool C_Object<dims>::operator!=(C_Object<dims>& obj) {
	return object_id != obj.object_id;
}

template <size_t dims>
inline bool C_Object<dims>::operator==(C_Object<dims>& obj) {
	return position == obj.position;
}

template <size_t dims>
inline void C_Object<dims>::reset() {
	buffer.reset();
}

template <size_t dims>
inline std::ostream& operator<<(std::ostream& output_stream, C_Object<dims>& obj) {
	output_stream << line
		<< "	Dane dla obiektu numer " << obj.object_id << '\n'
		<< "		Masa" << '\n'
		<< "			M = " << obj.mass << '\n'
		<< "		Promien" << '\n'
		<< "			R = " << obj.radius << '\n'
		<< "		Pozycja" << '\n'
		<< "			x = " << obj.position[0] << '\n'
		<< "			y = " << obj.position[1] << '\n'
		<< "			z = " << obj.position[2] << '\n'
		<< "		Predkosc" << '\n'
		<< "			x = " << obj.velocity[0] << '\n'
		<< "			y = " << obj.velocity[1] << '\n'
		<< "			z = " << obj.velocity[2] << '\n';
	return output_stream;
}

template <size_t dims>
inline std::istream& operator>>(std::istream& input_stream, C_Object<dims>& obj) {
	std::pair<double, int_fast32_t>
		fp_number = { double(),int_fast32_t() };
	auto
		to_double = [](std::pair<double, int_fast32_t> p)->double { return get<double>(p) * std::pow(10, get<int_fast32_t>(p)); };
	std::cout << line;
	std::cout << "Wprowadz parametry obiektu numer " << obj.object_id << "." << '\n';
	do {
		std::cout << "Podaj mase obiektu (kg): ";
		input_stream >> get<double>(fp_number) >> get<int_fast32_t>(fp_number);
	} while (input_error() || incorrect_value(get<double>(fp_number) < double()));
	obj.mass = to_double(fp_number);
	do {
		std::cout << "Podaj promien obiektu (m): ";
		input_stream >> get<double>(fp_number) >> get<int_fast32_t>(fp_number);
	} while (input_error() || incorrect_value(get<double>(fp_number) < double()));
	obj.radius = to_double(fp_number);
	do {
		input_error();
		std::cout << "Podaj polozenie (m): ";
		input_stream >> obj.position[0] >> obj.position[1] >> obj.position[2];
	} while (input_error());
	do {
		std::cout << "Podaj predkosc (m/s): ";
		input_stream >> obj.velocity[0] >> obj.velocity[1] >> obj.velocity[2];
	} while (input_error());
	return input_stream;
}

#endif // !C_obj