#include "C_main.h"

C_Main::C_Main() {
	t_counter.start();
	std::cout << line << "Space Physics Simulator - compiled: " << __DATE__ << " " << __TIME__ << '\n';
	event_log_obj() << "Program's initialization..." << _endl_;
	path_create();
	event_log_obj.create_log_file();
#ifdef _WIN32
	event_log_obj() << "Setting realtime priority class" << operation_evaluation<false>(Windows::SetPriorityClass(Windows::GetCurrentProcess(), ABOVE_NORMAL_PRIORITY_CLASS)) << _endl_;
#endif // _WIN32
	settings_obj = std::make_unique<C_Settings>(event_log_obj);
	event_log_obj() << "Initialization complete." << _endl_;
}

C_Main::~C_Main() {
	event_log_obj.finish_expected();
}

bool C_Main::run() {
	t_counter.stop();
	event_log_obj() << "Initialization time: " + std::to_string(t_counter.measured_timespan().count()) << _endl_;
	do {
		std::cout << line;
		std::cout << "Menu Glowne" << '\n';
		do {
			std::cout << line
				<< "0 - Zamknij program" << '\n'
				<< "1 - Symulacja" << '\n'
				<< "2 - Odczytywanie danych z analizy" << '\n'
				<< "3 - Ustawienia" << '\n'
				<< "4 - Obszar testowy" << '\n'
				<< "5 - Zresetuj program" << '\n'
				<< "Twoj wybor: ";
			std::cin >> main_choice;
		} while (input_error() || incorrect_value(main_choice > 5));
		event_log_obj() << "main_choice = " << main_choice << _endl_;
		switch (main_choice) {
		case 0:
			return false;
		case 1:
		{
			C_Universe<3>(settings_obj, event_log_obj).run();
			break;
		}
		case 2:
		{
			//running = std::make_unique<C_Analysis_Access>(settings_obj, event_log_obj);
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
	} while (std::clamp((int)main_choice, 1, 4) == main_choice);
	return true;
}

inline void C_Main::path_create() {
	if (!std::filesystem::exists(main_path)) {
		std::filesystem::create_directory(main_path);
		event_log_obj() << "Creating program's main path" << operation_evaluation<true>(std::filesystem::exists(main_path)) << _endl_;
	}
	std::filesystem::current_path(main_path);
	event_log_obj() << "Setting current directory" << operation_evaluation<true>(std::filesystem::current_path() == main_path) << _endl_;
	if (!std::filesystem::exists(main_path / L"simulations")) {
		std::filesystem::create_directory(main_path / L"simulations");
		event_log_obj() << "Creating program's simulation recording path" << operation_evaluation<true>(std::filesystem::exists(main_path / L"simulations")) << _endl_;
	}
	if (!std::filesystem::exists(main_path / L"logs")) {
		std::filesystem::create_directory(main_path / L"logs");
		event_log_obj() << "Creating program's event log path" << operation_evaluation<true>(std::filesystem::exists(main_path / L"logs")) << _endl_;;
	}
}