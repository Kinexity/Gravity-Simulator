#pragma once
#ifndef csdb
#define csbd
#include <cstdint>
#include <cstdlib>
#include <string>

struct S_Sim_Basic_Data {
	size_t
		num_of_objects = 0, //DO NOT CHANGE ORDER - 1
		sim_duration = 0, //DO NOT CHANGE ORDER - 2
		cycles_per_second_exponent = 0; //DO NOT CHANGE ORDER - 3
};

class C_Sim_Basic_Data {
private:
	S_Sim_Basic_Data
		internal_data;
public:
	explicit C_Sim_Basic_Data(S_Sim_Basic_Data sbd);
	C_Sim_Basic_Data() = default;
	~C_Sim_Basic_Data() = default;
	S_Sim_Basic_Data&
		operator()();
	void
		operator=(C_Sim_Basic_Data obj_to_assign);
	bool
		operator!=(C_Sim_Basic_Data obj_to_compare);
	static S_Sim_Basic_Data
		create_sbd(uint_fast64_t cycles_per_second_exponent, uint_fast64_t num_of_objects, uint_fast64_t sim_duration);
};

#endif // !csdb