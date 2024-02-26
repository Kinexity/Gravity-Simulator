#pragma once
#ifndef C_Mn
#define C_Mn
#include <iostream>
#include <string>
#include <fstream>
#include <limits>
#include <conio.h>
#include "C_Settings.h"
#include "C_Universe.h"
#include "C_Test.h"
#include "PureCPPLib/C_Event_Log.h"

class C_Main {
private:
	const std::wstring
		username = []()->auto{
		std::string temp = getenv("USERNAME");
		return std::wstring().assign(temp.begin(), temp.end());
	}();
	const std::filesystem::path
		main_path = L"C:/Users/" + username + L"/Documents/SpacePhysicsSimulator";
	uint_fast64_t
		main_choice = 0;
	std::unique_ptr<C_Settings>
		settings_obj;
	PCL::C_Event_Log
		event_log_obj{ true };
	PCL::C_Time_Counter
		t_counter;
public:
	C_Main();
	~C_Main();
	void
		path_create();
	bool
		run();
};

#endif // !C_Mn