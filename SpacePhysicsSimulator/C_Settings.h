#pragma once
#ifndef C_Set
#define C_Set
#include <fstream>
#include <filesystem>
#include <string>
#include <sstream>
#include "PureCPPLib/C_Event_Log.h"
#include "C_Sim_Basic_Data.h"
#include "constants.h"
using namespace std;
#undef max

struct S_Settings_Storage {
	uint_fast64_t
		num_of_sim = 0,
		num_of_threads_override = 0;
	double
		ram_GiB_limit = 1;
	bool
		enable_gpu_acceleration = bool(),
		full_cuda_acceleration = bool(),
		calculate_energy_error = bool();
	C_Sim_Basic_Data
		test_sim_basic_data = C_Sim_Basic_Data(C_Sim_Basic_Data::create_sbd(0, 4, 32));
};

class C_Settings {
private:
	mutex
		synchronizer;
	fstream
		settings_file;
	uint_fast64_t
		setting_choice = 0;
	S_Settings_Storage
		settings,
		settings_copy;
	PCL::C_Event_Log
		&event_log_obj;
	void
		read_settings(),
		update_settings(),
		log_settings();
public:
	void
		run();
	S_Settings_Storage&
		get_settings();
	template <typename Type, typename Func> void
		get_new_value(Type& variable, Func error_condition);
	C_Settings(PCL::C_Event_Log& ev_log_obj);
	~C_Settings();
};

template<typename Type, typename Func>
inline void C_Settings::get_new_value(Type & variable, Func error_condition) {
	do {
		std::cout << "Wartosc: ";
		cin >> variable;
	} while (input_error() || incorrect_value(error_condition()));
}

#endif // !C_Set