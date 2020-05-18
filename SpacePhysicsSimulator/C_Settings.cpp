#include "C_Settings.h"

#define change_check(x) if (settings_copy.x != settings.x) { event_log_obj() << "settings." << #x << ": " << settings_copy.x << " -> " << settings.x << _endl_; }

bool operator!=(S_Settings_Storage set_1, S_Settings_Storage set_2) {
	return
		(set_1.enable_gpu_acceleration != set_2.enable_gpu_acceleration) ||
		(set_1.num_of_sim != set_2.num_of_sim) ||
		(set_1.ram_GiB_limit != set_2.ram_GiB_limit) ||
		(set_1.num_of_threads_override != set_2.num_of_threads_override) ||
		(set_1.full_cuda_acceleration != set_2.full_cuda_acceleration) ||
		set_1.test_sim_basic_data != set_2.test_sim_basic_data;
}

C_Settings::C_Settings(PCL::C_Event_Log& ev_log_obj)
	: event_log_obj(ev_log_obj) {
	unique_lock<mutex> guard(synchronizer);
	settings_file.open("settings.bin", ios::out | ios::in | ios::binary);
	if (!settings_file.is_open()) {
		event_log_obj() << "Settings file doesn't exist! Creating new one..." << _endl_;
		settings_file.open("settings.bin", ios::out);
		settings_file.close();
		settings_file.open("settings.bin", ios::out | ios::in | ios::binary);
		event_log_obj() << "Creating new settings file..." << _endl_;
		settings_file.seekp(0, ios::beg);
		event_log_obj() << "Uploading default settings..." << _endl_;
		settings_file.write(reinterpret_cast<char*>(&settings), sizeof(settings));
	}
	event_log_obj() << "Opening settings file" << operation_evaluation<true>(settings_file.is_open() && settings_file) << _endl_;
	read_settings();
}

C_Settings::~C_Settings() {
	update_settings();
	settings_file.close();
	event_log_obj() << "Closing settings file" << operation_evaluation<false>(!settings_file.is_open()) << _endl_;
}

void C_Settings::read_settings() {
	settings_file.seekg(0, ios::beg);
	settings_file.read(reinterpret_cast<char*>(&settings), sizeof(settings));
	event_log_obj() << "Reading settings" << operation_evaluation<true>(settings_file.is_open() && settings_file) << _endl_;
	settings_copy = settings;
	log_settings();
}

void C_Settings::log_settings() {
	PCL::C_Event_Log_Buffer(event_log_obj, false, "SETTINGS_LOG")() << std::boolalpha << _endl_
		<< "Number of simulations: " << settings.num_of_sim << _endl_
		<< "Limit of the number of threads: " << settings.num_of_threads_override << _endl_
		<< "RAM usage limit (GiB): " << settings.ram_GiB_limit << _endl_
		<< "GPU acceleration enabled: " << settings.enable_gpu_acceleration << _endl_
		<< "Full CUDA acceleration enabled: " << settings.full_cuda_acceleration << _endl_ << std::noboolalpha;
}

void C_Settings::update_settings() {
	unique_lock<mutex> guard(synchronizer);
	if (settings != settings_copy) {
		change_check(num_of_sim);
		change_check(num_of_threads_override);
		change_check(ram_GiB_limit);
		change_check(enable_gpu_acceleration);
		change_check(full_cuda_acceleration);
		change_check(calculate_energy_error);
		change_check(test_sim_basic_data().num_of_objects);
		change_check(test_sim_basic_data().sim_duration);
		change_check(test_sim_basic_data().cycles_per_second_exponent);
		PCL::C_Event_Log_Buffer buffer(event_log_obj, true, "SETTINGS UPDATE");
		buffer() << "Updating settings..." << _endl_;
		settings_copy = settings;
		if (settings_file.is_open() && settings_file) {
			settings_file.seekp(0, ios::beg);
			settings_file.write(reinterpret_cast<char*>(&settings), sizeof(settings));
			if (settings_file) {
				buffer() << "Settings updated." << _endl_;
			}
			else {
				buffer() << "Error might have occured!" << _endl_;
				settings_file.clear();
			}
		}
		else {
			buffer() << "Settings file - unexpected error occured!" << _endl_;
		}
	}
}

S_Settings_Storage& C_Settings::get_settings() {
	unique_lock<mutex> guard(synchronizer);
	thread([&] {
		this_thread::sleep_for(500ms);
		update_settings();
	}).detach();
	return settings;
}

void C_Settings::run() {
	do {
		unique_lock<mutex> guard(synchronizer);
		std::cout << line << '\n';
		std::cout << "Menu Ustawien" << '\n';
		do {
			std::cout << line << '\n'
				<< "0 - Powrot do Glownego Menu" << '\n'
				<< "1 - Liczba symulacji (na wlasne ryzyko): " << settings.num_of_sim << '\n'
				<< "2 - Ograniczenie liczby watkow: " << settings.num_of_threads_override << '\n'
				<< "3 - Ograniczenie zuzycia RAMu: " << settings.ram_GiB_limit << '\n'
				<< "4 - Akceleracja na GPU: " << settings.enable_gpu_acceleration << '\n'
				<< "5 - Pelna akceleracja przy pomocy CUDA: " << settings.full_cuda_acceleration << '\n'
				<< "6 - Liczenie bledu energetycznego: " << settings.calculate_energy_error << '\n'
				<< "7 - Dane symulacji testowej:" << '\n'
				<< "	Liczba obiektow: " << settings.test_sim_basic_data().num_of_objects << '\n'
				<< "	Czas trwania: " << settings.test_sim_basic_data().sim_duration << '\n'
				<< "	Cykle na sekunde: " << settings.test_sim_basic_data().cycles_per_second_exponent << '\n'
				<< "Twoj wybor: ";
			cin >> setting_choice;
		} while (input_error() || incorrect_value(setting_choice > 7));
		event_log_obj() << "setting_choice = " << setting_choice << _endl_;
		if (setting_choice != 0) {
			std::cout << line << '\n';
		}
		{
			switch (setting_choice) {
			case 1: {
				get_new_value(settings.num_of_sim, [&]()->bool { return false; });
				break;
			}
			case 2: {
				do {
					std::cout << "Podaj liczbe watkow: ";
					cin >> settings.num_of_threads_override;
				} while (input_error());
				break;
			}
			case 3: {
				do {
					std::cout << "Podaj maksymalna ilosc RAMu (GiB): ";
					cin >> settings.ram_GiB_limit;
				} while (input_error());
				break;
			}
			case 4: {
				do {
					std::cout << "Akceleracja na GPU:" << '\n';
					std::cout << "0 - Blokuj" << '\n';
					std::cout << "1 - Uzywaj" << '\n';
					std::cout << "Twoj wybor: ";
					cin >> settings.enable_gpu_acceleration;
				} while (input_error());
				break;
			}
			case 5: {
				do {
					std::cout << "Pelna akceleracja przy pomocy CUDA: " << '\n';
					std::cout << "0 - Blokuj" << '\n';
					std::cout << "1 - Uzywaj" << '\n';
					std::cout << "Twoj wybor: ";
					cin >> settings.full_cuda_acceleration;
				} while (input_error());
				break;
			}
			case 6: {
				do {
					std::cout << "Liczenie bledu energetycznego: " << '\n';
					std::cout << "0 - Blokuj" << '\n';
					std::cout << "1 - Uzywaj" << '\n';
					std::cout << "Twoj wybor: ";
					cin >> settings.calculate_energy_error;
				} while (input_error());
				break;
			}
			case 7: {
				do {
					std::cout << "Podaj dane symulacji ([liczba obiektow : >= 2] [czas trwania : >= 1] [cykle na sekunde])" << '\n';
					std::cout << "Dane: ";
					cin >> settings.test_sim_basic_data().num_of_objects >> settings.test_sim_basic_data().sim_duration >> settings.test_sim_basic_data().cycles_per_second_exponent;
				} while (input_error() || incorrect_value(settings.test_sim_basic_data().num_of_objects < 2 || settings.test_sim_basic_data().sim_duration == 0));
			}
			}
			event_log_obj() << "User potentially changed settings" << _endl_;
		}
	} while (setting_choice != 0);
	update_settings();
}