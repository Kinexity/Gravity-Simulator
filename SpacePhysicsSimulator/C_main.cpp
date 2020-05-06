#include "C_main.h"

C_Main::C_Main() {
	t_counter.start();
	cout << line << endl << "Space Physics Simulator - Alpha 0.2" << endl;
	event_log_obj() << "Program's initialization..." << _endl_;
	path_create();
	event_log_obj.create_log_file();
#ifdef _WIN32
	event_log_obj() << "Setting realtime priority class" << operation_evaluation<false>(Windows::SetPriorityClass(Windows::GetCurrentProcess(), ABOVE_NORMAL_PRIORITY_CLASS)) << _endl_;
#endif // _WIN32
	settings_obj = make_unique<C_Settings>(event_log_obj);
	event_log_obj() << "Initialization complete." << _endl_;
}

C_Main::~C_Main() {
	event_log_obj.finish_expected();
}

bool C_Main::run() {
	t_counter.stop();
	event_log_obj() << "Initialization time: " + to_string(t_counter.measured_timespan().count()) << _endl_;
	do {
		cout << line << endl;
		cout << "Menu Glowne" << endl;
		do {
			cout << line << endl
				<< "0 - Zamknij program" << endl
				<< "1 - Symulacja" << endl
				<< "2 - Odczytywanie danych z analizy" << endl
				<< "3 - Ustawienia" << endl
				<< "4 - Obszar testowy" << endl
				<< "5 - Zresetuj program" << endl
				<< "Twoj wybor: ";
			cin >> main_choice;
		} while (input_error() || incorrect_value(main_choice > 5));
		event_log_obj() << "main_choice = " << main_choice << _endl_;
		switch (main_choice) {
		case 0:
			return false;
		case 1:
		{
			C_Universe(settings_obj, event_log_obj).run();
			break;
		}
		case 2:
		{
			//running = make_unique<C_Analysis_Access>(settings_obj, event_log_obj);
			break;
		}
		case 3:
		{
			settings_obj->run();
			break;
		}
		case 4:
		{
			C_Test(event_log_obj).run();
			break;
		}
		}
		if (running != nullptr) {
			running->run();
			running.reset();
		}
	} while (PCL::math::interval<decltype(main_choice), true, true>(1,4) == main_choice);
	return true;
}

inline void C_Main::path_create() {
	if (!filesystem::exists(main_path)) {
		filesystem::create_directory(main_path);
		event_log_obj() << "Creating program's main path" << operation_evaluation<true>(filesystem::exists(main_path)) << _endl_;
	}
	filesystem::current_path(main_path);
	bool b = (filesystem::current_path() == main_path);
	event_log_obj() << "Setting current directory" << operation_evaluation<true>(b) << _endl_;
	if (!filesystem::exists(main_path / L"simulations")) {
		filesystem::create_directory(main_path / L"simulations");
		event_log_obj() << "Creating program's simulation recording path" << operation_evaluation<true>(filesystem::exists(main_path / L"simulations")) << _endl_;
	}
	if (!filesystem::exists(main_path / L"logs")) {
		filesystem::create_directory(main_path / L"logs");
		event_log_obj() << "Creating program's event log path" << operation_evaluation<true>(filesystem::exists(main_path / L"logs")) << _endl_;;
	}
}