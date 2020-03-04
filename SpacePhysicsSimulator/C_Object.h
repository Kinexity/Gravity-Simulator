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
#include "C_Universe.h"
#include "C_Vector.h"
#include "C_AVX_Vector.h"
using namespace std;

class C_Object {
private:
	//C_Universe&
		//parent_universe;
	std::array<double, 3>
		position,
		velocity,
		acceleration = decltype(C_Object::position){0};
	double
		mass = double(),
		radius = double(),
		standard_grav_param = double();
	unique_ptr<decltype(C_Object::position)[]>
		buffer;
public:
	uint_fast64_t
		object_id = 0;
	//C_Object(C_Universe& universe_ref);
	C_Object() = default;
	C_Object(C_Object& obj);
	C_Object(C_Object&&) = default;
	C_Object& operator=(const C_Object&) = default;
	~C_Object() = default;
	bool
		operator!=(C_Object& obj),
		operator==(C_Object& obj);
	void
		update(uint_fast64_t sec_ord),
		reset(),
		allocate(uint_fast64_t period_lenght),
		data_rand(uint_fast64_t noo),
		calculate_sgp();
	double
		operator()(C_Object& obj_2),
		velocity_val();
	friend class
		C_Universe;
	friend class
		C_Analysis;
	friend ostream& operator<<(ostream& str, C_Object& obj);
	friend istream& operator>>(istream& str, C_Object& obj);
};

#endif // !C_obj