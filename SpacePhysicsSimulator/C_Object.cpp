#include "C_Object.h"

void C_Object::calculate_sgp() {
	standard_grav_param = mass * 6.67408313131313131313e-11;
}

double C_Object::operator()(C_Object& obj_2) {
	return std::hypot(position[0] - obj_2.position[0], position[1] - obj_2.position[1], position[2] - obj_2.position[2]);
}

double C_Object::velocity_val() {
	return std::hypot(velocity[0], velocity[1], velocity[2]);
}

void C_Object::data_rand(uint_fast64_t noo) {
	srand(clock() * object_id);
	mass = static_cast<double>(rnd() % 900 + 100)* pow(10, rnd() % 20 + 10);
	radius = pow(10, 1);
	auto x = static_cast<uint_fast64_t>(sqrt(cbrt(noo)));
	auto min_pow = 6;
	position = {
		((static_cast<double>(rnd() % 9000000 + 1000000))* pow(10, rnd() % 5 + min_pow)* ((rnd() % 2 == 1) ? 1.0 : -1.0)),
		((static_cast<double>(rnd() % 9000000 + 1000000))* pow(10, rnd() % 5 + min_pow)* ((rnd() % 2 == 1) ? 1.0 : -1.0)),
		((static_cast<double>(rnd() % 9000000 + 1000000))* pow(10, rnd() % 5 + min_pow)* ((rnd() % 2 == 1) ? 1.0 : -1.0))
	};
	velocity = {
		(static_cast<double>(rnd() % 9000 + 1000)* ((rnd() % 2 == 1) ? 1.0 : -1.0)),
		(static_cast<double>(rnd() % 9000 + 1000)* ((rnd() % 2 == 1) ? 1.0 : -1.0)),
		(static_cast<double>(rnd() % 9000 + 1000)* ((rnd() % 2 == 1) ? 1.0 : -1.0))
	};
}

//C_Object::C_Object(C_Universe& universe_ref) : 
//parent_universe(universe_ref) {}

C_Object::C_Object(C_Object& obj) {
	memcpy(&position, &obj.position, ((int)&buffer - (int)&position));
}

bool C_Object::operator!=(C_Object& obj) {
	return object_id != obj.object_id;
}

bool C_Object::operator==(C_Object& obj) {
	return position == obj.position;
}

void C_Object::update(uint_fast64_t sec_ord) {
	buffer[sec_ord] = position;
}

void C_Object::allocate(uint_fast64_t period_lenght) {
	buffer = make_unique<decltype(C_Object::position)[]>(period_lenght + 1);
}

void C_Object::reset() {
	buffer.reset();
}

ostream& operator<<(ostream& output_stream, C_Object& obj) {
	output_stream << line << endl
		<< "	Dane dla obiektu numer " << obj.object_id << endl
		<< "		Masa" << endl
		<< "			M = " << obj.mass << endl
		<< "		Promien" << endl
		<< "			R = " << obj.radius << endl
		<< "		Pozycja" << endl
		<< "			x = " << obj.position[0] << endl
		<< "			y = " << obj.position[1] << endl
		<< "			z = " << obj.position[2] << endl
		<< "		Predkosc" << endl
		<< "			x = " << obj.velocity[0] << endl
		<< "			y = " << obj.velocity[1] << endl
		<< "			z = " << obj.velocity[2] << endl;
	return output_stream;
}

istream& operator>>(istream& input_stream, C_Object& obj) {
	pair<double, int_fast32_t>
		fp_number = { double(),int_fast32_t() };
	auto
		to_double = [](pair<double, int_fast32_t> p)->double { return get<double>(p) * pow(10, get<int_fast32_t>(p)); };
	cout << line << endl;
	cout << "Wprowadz parametry obiektu numer " << obj.object_id << "." << endl;
	do {
		cout << "Podaj mase obiektu (kg): ";
		input_stream >> get<double>(fp_number) >> get<int_fast32_t>(fp_number);
	} while (input_error() || incorrect_value(get<double>(fp_number) < double()));
	obj.mass = to_double(fp_number);
	do {
		cout << "Podaj promien obiektu (m): ";
		input_stream >> get<double>(fp_number) >> get<int_fast32_t>(fp_number);
	} while (input_error() || incorrect_value(get<double>(fp_number) < double()));
	obj.radius = to_double(fp_number);
	do {
		input_error();
		cout << "Podaj polozenie (m): ";
		input_stream >> obj.position[0] >> obj.position[1] >> obj.position[2];
	} while (input_error());
	do {
		cout << "Podaj predkosc (m/s): ";
		input_stream >> obj.velocity[0] >> obj.velocity[1] >> obj.velocity[2];
	} while (input_error());
	return input_stream;
}