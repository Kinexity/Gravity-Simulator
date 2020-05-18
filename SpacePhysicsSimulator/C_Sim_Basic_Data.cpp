#include "C_Sim_Basic_Data.h"

C_Sim_Basic_Data::C_Sim_Basic_Data(S_Sim_Basic_Data sbd) : internal_data(sbd) {}

S_Sim_Basic_Data& C_Sim_Basic_Data::operator()() {
	return internal_data;
}

void C_Sim_Basic_Data::operator=(C_Sim_Basic_Data obj_to_assign) {
	std::memcpy(this, &obj_to_assign, sizeof(*this));
}

bool C_Sim_Basic_Data::operator!=(C_Sim_Basic_Data obj_to_compare) {
	return internal_data.num_of_objects != obj_to_compare.internal_data.num_of_objects ||
		internal_data.sim_duration != obj_to_compare.internal_data.sim_duration ||
		internal_data.cycles_per_second_exponent != obj_to_compare.internal_data.cycles_per_second_exponent;
}

S_Sim_Basic_Data C_Sim_Basic_Data::create_sbd(uint_fast64_t cycles_per_second_exponent, uint_fast64_t num_of_objects, uint_fast64_t sim_duration) {
	S_Sim_Basic_Data s;
	s.cycles_per_second_exponent = cycles_per_second_exponent;
	s.num_of_objects = num_of_objects;
	s.sim_duration = sim_duration;
	return s;
}