#pragma once
#define _ENABLE_EXTENDED_ALIGNED_STORAGE
#ifndef C_univ
#define C_univ
#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <execution>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <thread>
#include <vector>
#include <immintrin.h>
#include <functional>
#include <set>
#include <map>
#include <list>
#include "PureCPPLib/C_Barrier.h"
#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN 
namespace Windows {
#include <Windows.h>
}
#endif // _WIN32
#include "C_Settings.h"
#include "PureCPPLib/C_thread_set.h"
#include "C_Base.h"
#include "C_Object.h"
#include "PureCPPLib/C_Event_Log.h"
#include "PureCPPLib/C_Indexer.h"
#include "PureCPPLib/C_Time_Counter.h"
#include "PureCPPLib/C_Random.h"
#include "constants.h"
#undef max
#undef min

class C_Settings;

inline constexpr auto dimensions = 3;

class C_Universe : public C_Base {
private:
	uint_fast64_t
		sim_ord = 0,
		action_identifier = 0,
		num_of_threads_in_use = 0,
		range_size = 0,
		period_lenght = 0,
		periods = 0,
		mod_choice = 0;
	C_Sim_Basic_Data
		sim_basic_data;
	fstream
		backup_file;
	mutex
		mut_ex;
	const std::filesystem::path
		recording_path = filesystem::current_path() / L"simulations";
	unique_ptr<PCL::C_Barrier>
		barrier;
	std::vector<C_Object>
		obj_arr;
	unique_ptr<C_Settings>&
		settings_obj;
	PCL::C_Event_Log&
		event_log_obj; 
	unique_ptr<PCL::C_Event_Log_Buffer>
		event_log_buffer;
	unique_ptr<PCL::C_Indexer<>>
		indexer;
	PCL::C_Time_Counter
		time_counter;
	bool
		bool_pre_sim_edit = false,
		abort_simulation = false,
		include_analysis = false;
	std::once_flag
		sync_flags[15];
	std::atomic<double>
		E_c_0 = 0,
		E_c_k = 0;
public:
	C_Universe(unique_ptr<C_Settings>& set_obj, PCL::C_Event_Log& ev_log_obj);
	~C_Universe() = default;
	void
		run() override;
	inline void
		data_read(),
		data_write(),
		data_main_input(),
		data_speed_test(),
		data_display(),
		data_change(),
		data_log(),
		simulation_directory_create(),
		simulation_packed_thread(), 
		simulation_packed(),
		simulation_preperator(),
		action_choice(),
		result_recording();
};

#endif // !C_univ