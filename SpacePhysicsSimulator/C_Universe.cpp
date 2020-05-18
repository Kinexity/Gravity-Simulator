#include "C_Universe.h"
using namespace placeholders;

inline void C_Universe::action_choice() {
	do {
		std::cout << line << '\n';
		std::cout << "Wybierz sposob dzialania: " << '\n';
		std::cout << "0 - Powrot do Menu" << '\n';
		std::cout << "1 - Podaj dane i uruchom symulacje" << '\n';
		std::cout << "2 - Wczytaj z pliku" << '\n';
		std::cout << "3 - Utworz plik do wczytania" << '\n';
		std::cout << "4 - Test wydajnosci" << '\n';
		std::cout << "Twoj wybor: ";
		cin >> action_identifier;
	} while (input_error() || incorrect_value(action_identifier > 4));
	event_log_obj() << "action_identifier = " << action_identifier << _endl_;
}

C_Universe::C_Universe(unique_ptr<C_Settings>& set_obj, PCL::C_Event_Log& ev_log_obj) :
	settings_obj(set_obj),
	event_log_obj(ev_log_obj) {
	event_log_obj() << "Simulation write path: " << recording_path.string() << _endl_;
	action_choice();
	PCL::C_Indexer object_indexer(sim_basic_data().num_of_objects);
	switch (action_identifier) {
	case 0:
	{
		abort_simulation = true;
		break;
	}
	case 1:
	{
		data_main_input();
		obj_arr.resize(sim_basic_data().num_of_objects);
		for (C_Object<dimensions>& object : obj_arr) {
			object.object_id = object_indexer();
			//std::cin >> object;
		}
		break;
	}
	case 2:
	{
		data_read();
		if (!abort_simulation) {
			for (C_Object<dimensions>& object : obj_arr) {
				object.object_id = object_indexer();
			}
		}
		break;
	}
	case 3:
	{
		data_main_input();
		obj_arr.resize(sim_basic_data().num_of_objects);
		for (C_Object<dimensions>& object : obj_arr) {
			object.object_id = object_indexer();
			//cin >> object;
		}
		abort_simulation = true;
		break;
	}
	case 4:
	{
		std::set<decltype(C_Object<dimensions>::position)> pos_set;
		uint_fast64_t duplicates_count = 0;
		data_speed_test();
		obj_arr.resize(sim_basic_data().num_of_objects);
		for (auto& object : obj_arr) {
			object.object_id = object_indexer();
			do {
				duplicates_count++;
				object.data_rand(sim_basic_data().num_of_objects);
				pos_set.insert(object.position);
			} while (pos_set.size() != object.object_id + 1);
		}
		duplicates_count -= sim_basic_data().num_of_objects;
		event_log_obj() << "Object random gen position duplicates count: " << duplicates_count << _endl_;
	}
	}
	if (action_identifier != 0) {
		event_log_obj() << "Objects' array initialized." << _endl_;
		if (action_identifier != 3) {
			period_lenght = std::min(static_cast<uint_fast64_t>(
				pow(
					2,
					static_cast<uint_fast64_t>(
						floor(
							log2(
								settings_obj->get_settings().ram_GiB_limit * pow(2, 30) / sizeof(decltype(C_Object<dimensions>::position)) / sim_basic_data().num_of_objects
							)
						)
						)
				)
				), sim_basic_data().sim_duration);
			if (period_lenght > 0) {
				if (period_lenght != sim_basic_data().sim_duration) {
					periods = static_cast<uint_fast64_t>(ceil(static_cast<double>(sim_basic_data().sim_duration) / period_lenght));
				}
				else {
					periods = 1;
				}
				data_log();
				event_log_obj() << "Necessery memory acquired." << _endl_;
			}
			else {
				event_log_obj() << "Not enough RAM! Aborting simulation!";
				abort_simulation = true;
			}
		}
	}
}

void C_Universe::run() {
	if (!abort_simulation) {
		event_log_obj() << (action_identifier != 4 ? "Simulation" : "Speed test") << " started." << _endl_;
		simulation_directory_create();
		data_change();
		data_write();
		simulation_preperator();
		std::cout << fixed
			<< line << '\n'
			<< "Czas trwania symulacji: " << time_counter.measured_timespan().count() << '\n'
			<< line << '\n'
			<< "Czas trwania symulacji jednego okresu: " << time_counter.measured_timespan().count() / periods << '\n'
			<< scientific
			<< "Oddzialywania na sekunde: " << ((sim_basic_data().num_of_objects * (sim_basic_data().num_of_objects - 1) * sim_basic_data().sim_duration) / time_counter.measured_timespan().count() / 2) << '\n'
			<< line << '\n';
		if (settings_obj->get_settings().calculate_energy_error) {
			std::cout << "Blad energetyczny: " << (abs(E_c_k / E_c_0) - 1) << '\n';
		}
		event_log_obj()
			<< (action_identifier != 4 ? "Simulation" : "Speed test") << " finished - (1 period / all periods) " << time_counter.measured_timespan().count() / sim_basic_data().sim_duration << " / " << time_counter.measured_timespan().count() << _endl_
			<< "Gravitational interactions per second: " << ((sim_basic_data().num_of_objects * (sim_basic_data().num_of_objects - 1) * sim_basic_data().sim_duration) / time_counter.measured_timespan().count() / 2) << _endl_;
		if (settings_obj->get_settings().calculate_energy_error) {
			event_log_obj()
				<< "Energy_0 = " << (double)E_c_0 << _endl_
				<< "Energy_k = " << (double)E_c_k << _endl_
				<< "Energy error: " << fixed << (abs((double)E_c_k / (double)E_c_0) - 1) << scientific << " %" << _endl_;
		}
		result_recording();
		settings_obj->get_settings().num_of_sim++;
	}
	else {
		std::cout << line << '\n' << "Symulacja " << (action_identifier == 0 || action_identifier == 3 ? "pominieta" : "anulowana") << '\n';
		event_log_obj() << "Simulation " << (action_identifier == 0 || action_identifier == 3 ? "cancelled" : "aborted") << _endl_;
	}
}

inline void C_Universe::simulation_directory_create() {
	auto simulation_directory_path = recording_path / (std::wstring(L"/simulation[") + to_wstring(settings_obj->get_settings().num_of_sim) + L"]" + std::wstring(action_identifier == 4 ? L"[ST]" : L""));
	event_log_obj() << "Creating simulation's directory" << operation_evaluation<false>(std::filesystem::create_directory(simulation_directory_path)) << _endl_;
}

inline void C_Universe::data_write() {
	backup_file.open(recording_path / (std::wstring(L"simulation[") + to_wstring(settings_obj->get_settings().num_of_sim) + L"]" + std::wstring(action_identifier == 4 ? L"[ST]" : L"")) / L"objects.bin", ios::out | ios::app | ios::binary);
	backup_file.write(reinterpret_cast<char*>(&sim_basic_data()), sizeof(sim_basic_data()));
	for_each(execution::seq, obj_arr.begin(), obj_arr.end(), [&](C_Object<dimensions>& obj)->void {
		backup_file.write(reinterpret_cast<char*>(&obj.position), (sizeof(decltype(obj.position))));
		backup_file.write(reinterpret_cast<char*>(&obj.velocity), (sizeof(decltype(obj.velocity))));
		backup_file.write(reinterpret_cast<char*>(&obj.mass), (sizeof(decltype(obj.mass))));
		backup_file.write(reinterpret_cast<char*>(&obj.radius), (sizeof(decltype(obj.radius))));
	});
	event_log_obj() << "Saving simulation's details to file" << operation_evaluation<false>(backup_file.is_open()) << _endl_;
	backup_file.close();
}

inline void C_Universe::data_read() {
	do {
		std::cout << line << '\n';
		std::cout << "Liczba symulacji: " << settings_obj->get_settings().num_of_sim << '\n';
		std::cout << line << '\n';
		std::cout << "Podaj numer symulacji: ";
		cin >> sim_ord;
	} while (input_error() || incorrect_value(sim_ord >= settings_obj->get_settings().num_of_sim));
	if (backup_file.is_open()) {
		backup_file.close();
	}
	backup_file.open("simulations/simulation[" + to_string(sim_ord) + "]/objects.bin", ios::in | ios::binary);
	if (backup_file.is_open()) {
		backup_file.seekg(0, ios::beg); // zmiana
		backup_file.read(reinterpret_cast<char*>(&sim_basic_data()), sizeof(sim_basic_data()));
		obj_arr.resize(sim_basic_data().num_of_objects);
		for_each(execution::par, obj_arr.begin(), obj_arr.end(), [&](C_Object<dimensions>& obj)->void {
			backup_file.read(reinterpret_cast<char*>(&obj.position), (sizeof(decltype(C_Object<dimensions>::position))));
			backup_file.read(reinterpret_cast<char*>(&obj.velocity), (sizeof(decltype(C_Object<dimensions>::velocity))));
			backup_file.read(reinterpret_cast<char*>(&obj.mass), (sizeof(decltype(C_Object<dimensions>::mass))));
			backup_file.read(reinterpret_cast<char*>(&obj.radius), (sizeof(decltype(C_Object<dimensions>::radius))));
		});
		event_log_obj() << "Reading simulation's details from file" << operation_evaluation<false>(static_cast<bool>(backup_file)) << _endl_;
		abort_simulation |= !backup_file;
		backup_file.close();
		data_display();
	}
	else {
		event_log_obj() << "Couldn't open file: " << filesystem::current_path() / (string("simulations/simulation[") + to_string(sim_ord) + "]/objects.bin") << _endl_;
		abort_simulation = true;
	}
}

inline void C_Universe::simulation_packed_thread() {
	uint32_t
		sync_flag_index = 0;
	static unique_ptr< unique_ptr<double[]>[]>
		position,
		velocity,
		acceleration;
	static unique_ptr<double[]>
		sgp_arr;
	static unique_ptr<mutex[]>
		mutex_arr;
	const auto
		thr_index = (*indexer)();
	barrier->arrive_and_wait();
	const double
		work = static_cast<double>(sim_basic_data().num_of_objects) * (sim_basic_data().num_of_objects - 1) / (2 * num_of_threads_in_use);
	auto obj_arr_index = [&](uint_fast64_t thr_index) {
		uint_fast64_t
			ret_index = 0;
		double
			work_to_index = double(),
			work_st = thr_index * work,
			min_diff = numeric_limits<double>::infinity();
		if (thr_index == num_of_threads_in_use) {
			return sim_basic_data().num_of_objects;
		}
		for (uint_fast64_t index = 0; index < sim_basic_data().num_of_objects; index += 4) {
			if (abs(work_st - work_to_index) < min_diff) {
				ret_index = index;
				min_diff = abs(work_st - work_to_index);
			}
			else {
				return ret_index;
			}
			work_to_index += 2 * ((sim_basic_data().num_of_objects - index - 1) + (sim_basic_data().num_of_objects - index - 4));
		}
		return ret_index;
	};
	auto to_4_mul = [](uint_fast64_t num) { return num - num % 4 + (num % 4 > 0 ? 4 : 0); };
	const uint_fast64_t
		arr_start_point = obj_arr_index(thr_index),
		arr_end_point = obj_arr_index(thr_index + 1),
		eq_arr_start_point = to_4_mul(sim_basic_data().num_of_objects * thr_index / num_of_threads_in_use),
		eq_arr_end_point = to_4_mul(sim_basic_data().num_of_objects * (thr_index + 1) / num_of_threads_in_use);
	static std::vector<std::reference_wrapper<decltype(arr_start_point)>> arr_start_point_arr;
	{
		static std::list<std::reference_wrapper<decltype(arr_start_point)>> arr_start_point_arr_temp;
		std::call_once(sync_flags[sync_flag_index++], [&] {
			arr_start_point_arr_temp.clear();
		});
		barrier->arrive_and_wait();
		{
			std::unique_lock<std::mutex> lck(mut_ex);
			arr_start_point_arr_temp.push_back(arr_start_point);
		}
		barrier->arrive_and_wait();
		std::call_once(sync_flags[sync_flag_index++], [&] {
			arr_start_point_arr = { arr_start_point_arr_temp.begin(), arr_start_point_arr_temp.end() };
			std::sort(arr_start_point_arr.begin(), arr_start_point_arr.end());
		});
		barrier->arrive_and_wait();
	}
	std::call_once(sync_flags[sync_flag_index++], [&] {
		position = std::make_unique<unique_ptr<double[]>[]>(dimensions);
		velocity = std::make_unique<unique_ptr<double[]>[]>(dimensions);
		acceleration = std::make_unique<unique_ptr<double[]>[]>(dimensions);
		for (auto axis = 0; axis < dimensions; axis++) {
			position[axis] = std::make_unique<double[]>(sim_basic_data().num_of_objects + simd_width - 1);
			velocity[axis] = std::make_unique<double[]>(sim_basic_data().num_of_objects + simd_width - 1);
			acceleration[axis] = std::make_unique<double[]>(sim_basic_data().num_of_objects + simd_width - 1);
			std::for_each(std::execution::par_unseq, obj_arr.begin(), obj_arr.end(), [&](C_Object<dimensions>& obj) {
				position[axis][obj.object_id] = obj.position[axis];
				velocity[axis][obj.object_id] = obj.velocity[axis];
			});
		}
		sgp_arr = std::make_unique<double[]>(sim_basic_data().num_of_objects + simd_width - 1);
		std::for_each(std::execution::par_unseq, obj_arr.begin(), obj_arr.end(), [&](C_Object<dimensions>& obj) {
			sgp_arr[obj.object_id] = obj.standard_grav_param;
		});
		for (auto obj_id = sim_basic_data().num_of_objects; obj_id < sim_basic_data().num_of_objects + simd_width - 1; obj_id++) {
			sgp_arr[obj_id] = 0;
			for (auto axis = 0; axis < dimensions; axis++) {
				position[axis][obj_id] = std::pow(rnd() % 10, rnd() % 10 + 10);
			}
		}
		mutex_arr = std::make_unique<mutex[]>(sim_basic_data().num_of_objects + simd_width - 1);
	});
	barrier->arrive_and_wait();
	{ //pre-sim log
		std::call_once(sync_flags[sync_flag_index++], [&] {
			event_log_buffer = std::make_unique<PCL::C_Event_Log_Buffer>(event_log_obj, true, "THREAD LOG");
		});
		barrier->arrive_and_wait();
		{
			std::unique_lock<std::mutex> lck(mut_ex);
			*event_log_buffer << "Thread " << thr_index << " - <" << arr_start_point << ";" << arr_end_point << ")" << _endl_;
		}
		barrier->arrive_and_wait();
		std::call_once(sync_flags[sync_flag_index++], [&] {
			event_log_buffer = std::make_unique<PCL::C_Event_Log_Buffer>(event_log_obj, true, "SIMULATION TIME LOG");
			*event_log_buffer << "Period	G_Time	A_Time" << _endl_;
		});
	} //pre-sim log end
	auto cycles_per_second = std::pow(2, sim_basic_data().cycles_per_second_exponent);
	auto cps = _mm256_set1_pd(double(1) / cycles_per_second);
	auto lock = [&](uint_fast64_t index) {
		mutex_arr[index].lock();
		mutex_arr[index + 1].lock();
		mutex_arr[index + 2].lock();
		mutex_arr[index + 3].lock();
	};
	auto unlock = [&](uint_fast64_t index) {
		mutex_arr[index].unlock();
		mutex_arr[index + 1].unlock();
		mutex_arr[index + 2].unlock();
		mutex_arr[index + 3].unlock();
	};
	auto system_energy = [&] {
		return std::transform_reduce(std::execution::par_unseq, obj_arr.begin() + thr_index * obj_arr.size() / num_of_threads_in_use, obj_arr.begin() + (thr_index + 1) * obj_arr.size() / num_of_threads_in_use, 0.0, std::plus<>(), [&](auto& obj) {
			return obj.mass *
				(obj.velocity_val() / 2 -
					std::transform_reduce(std::execution::par_unseq, obj_arr.begin() + obj.object_id + 1, obj_arr.end(), 0.0, std::plus<>(), [&](auto& op_obj) {
				return op_obj.standard_grav_param / op_obj(obj);
			}));
		}
		);
	};
	const auto
		zero_avx = _mm256_setzero_pd();
	__m256d
		pos[dimensions] = { zero_avx },
		vec[dimensions] = { zero_avx },
		acc_obj[dimensions] = { zero_avx },
		sgp = zero_avx;
	auto
		temp_acceleration = std::make_unique<unique_ptr<double[]>[]>(dimensions);
	barrier->arrive_and_wait();
	static std::vector<uint_fast64_t> range_lims;
	std::call_once(sync_flags[sync_flag_index++], [&] {
		auto work = arr_start_point_arr.size() * sim_basic_data().num_of_objects - std::reduce(arr_start_point_arr.begin(), arr_start_point_arr.end(), uint_fast64_t());
		work = work / num_of_threads_in_use + (work % num_of_threads_in_use == 0 ? 0 : 1);
		range_lims.resize(num_of_threads_in_use + 1);
		range_lims[0] = 0;
		size_t
			lim_index = 1,
			obj_it_index = 0,
			work_count = 0,
			thr_work_collect = 1;
		for (size_t it_thr_index = 1; it_thr_index < num_of_threads_in_use; it_thr_index++) {
			bool calculated = true;
			do {

				obj_it_index += 4;
				if (obj_it_index >= arr_start_point_arr[lim_index]) {
					//thr_work_collect
				}
			} while (!calculated);
		}
		range_lims[num_of_threads_in_use] = sim_basic_data().num_of_objects;
	});
	for (auto axis = 0; axis < dimensions; axis++) {
		temp_acceleration[axis] = std::make_unique<double[]>(sim_basic_data().num_of_objects + simd_width - 1);
	}
	if (settings_obj->get_settings().calculate_energy_error) {
		barrier->arrive_and_wait();
		{
			auto temp_energy = system_energy();
			std::unique_lock<std::mutex> lck(mut_ex);
			E_c_0 += temp_energy;
		}
	}
	barrier->arrive_and_wait();
	std::call_once(sync_flags[sync_flag_index++], [&] {
		if (settings_obj->get_settings().calculate_energy_error) {
			event_log_obj() << nm(E_c_0) << _endl_;
		}
		time_counter.start();
	});
	barrier->arrive_and_wait();
	for (uint_fast64_t period_ord = 0; period_ord < periods; period_ord++) {
		if (thr_index == 0) {
			tc.start();
			if (period_ord == periods - 1) {
				period_lenght = sim_basic_data().sim_duration - (periods - 1) * period_lenght;
			}
		}
		barrier->arrive_and_wait();
		for (uint_fast64_t second = 1; second <= period_lenght; second++) {
			for (uint_fast64_t sec_cycle = 0; sec_cycle < cycles_per_second; sec_cycle++) {
				for (auto axis = 0; axis < dimensions; axis++) {
					std::fill(std::execution::par, &(temp_acceleration[axis][arr_start_point]), &(temp_acceleration[axis][sim_basic_data().num_of_objects]), 0.0); //clear temp_acceleration
					std::fill(std::execution::par, &(acceleration[axis][eq_arr_start_point]), &(acceleration[axis][eq_arr_end_point]), 0.0); //clear acceleration
				}
				barrier->arrive_and_wait();
				for (uint_fast64_t obj_id = arr_start_point; obj_id < arr_end_point; obj_id += 4) {
					sgp = _mm256_load_pd(&(sgp_arr[obj_id]));
					for (auto axis = 0; axis < dimensions; axis++) {
						pos[axis] = _mm256_load_pd(&(position[axis][obj_id]));
						acc_obj[axis] = zero_avx;
					}
					for (uint_fast64_t op_obj_id = obj_id + 1; op_obj_id < arr_end_point; op_obj_id++) { //objects assigned to this thread
						auto dist_2 = zero_avx;
						for (auto axis = 0; axis < dimensions; axis++) {
							vec[axis] = _mm256_sub_pd(_mm256_load_pd(&(position[axis][op_obj_id])), pos[axis]); //calculate distance vectors
							dist_2 = _mm256_fmadd_pd(vec[axis], vec[axis], dist_2);
						}
						const auto&& dist_3 = _mm256_mul_pd(dist_2, _mm256_sqrt_pd(dist_2));
						const auto&& scaled_sgp_op = _mm256_div_pd(_mm256_load_pd(&(sgp_arr[op_obj_id])), dist_3);
						const auto&& scaled_sgp = _mm256_div_pd(sgp, dist_3);
						for (auto axis = 0; axis < dimensions; axis++) {
							acc_obj[axis] = _mm256_fmadd_pd(vec[axis], scaled_sgp_op, acc_obj[axis]); //obj_id acceleration from op_obj_id
							_mm256_store_pd(&(acceleration[axis][op_obj_id]), _mm256_fnmadd_pd(vec[axis], scaled_sgp, _mm256_load_pd(&(acceleration[axis][op_obj_id])))); //op_obj_id acceleration from obj_id and save to temp_acceleration
						}
					}
					for (uint_fast64_t op_obj_id = arr_end_point; op_obj_id < sim_basic_data().num_of_objects; op_obj_id++) { //objects NOT assigned to this thread
						auto dist_2 = zero_avx;
						for (auto axis = 0; axis < dimensions; axis++) {
							vec[axis] = _mm256_sub_pd(_mm256_load_pd(&(position[axis][op_obj_id])), pos[axis]); //calculate distance vectors
							dist_2 = _mm256_fmadd_pd(vec[axis], vec[axis], dist_2);
						}
						const auto&& dist_3 = _mm256_mul_pd(dist_2, _mm256_sqrt_pd(dist_2));
						const auto&& scaled_sgp_op = _mm256_div_pd(_mm256_load_pd(&(sgp_arr[op_obj_id])), dist_3);
						const auto&& scaled_sgp = _mm256_div_pd(sgp, dist_3);
						for (auto axis = 0; axis < dimensions; axis++) {
							acc_obj[axis] = _mm256_fmadd_pd(vec[axis], scaled_sgp_op, acc_obj[axis]); //obj_id acceleration from op_obj_id
							_mm256_store_pd(&(temp_acceleration[axis][op_obj_id]), _mm256_fnmadd_pd(vec[axis], scaled_sgp, _mm256_load_pd(&(temp_acceleration[axis][op_obj_id])))); //op_obj_id acceleration from obj_id and save to temp_acceleration
						}
					}
					for (auto axis = 0; axis < dimensions; axis++) {
						_mm256_store_pd(&(acceleration[axis][obj_id]), _mm256_add_pd(acc_obj[axis], _mm256_load_pd(&(acceleration[axis][obj_id])))); //obj_id acceleration save to temp_acceleration
					}
				}
				barrier->arrive_and_wait();
				for (auto axis = 0; axis < dimensions; axis++) {
					for (uint_fast64_t op_obj_id = arr_end_point; op_obj_id < sim_basic_data().num_of_objects; op_obj_id += 4) {
						lock(op_obj_id);
						_mm256_store_pd(&(acceleration[axis][op_obj_id]), _mm256_add_pd(_mm256_load_pd(&(acceleration[axis][op_obj_id])), _mm256_load_pd(&(temp_acceleration[axis][op_obj_id])))); //sum acceleration from threads
						unlock(op_obj_id);
					}
				}
				barrier->arrive_and_wait();
				for (auto axis = 0; axis < dimensions; axis++) {
					for (auto obj_id = eq_arr_start_point; obj_id < eq_arr_end_point; obj_id += 4) {
						_mm256_store_pd(&(velocity[axis][obj_id]), _mm256_fmadd_pd(_mm256_load_pd(&(acceleration[axis][obj_id])), cps, _mm256_load_pd(&(velocity[axis][obj_id])))); //update position and velocity
						_mm256_store_pd(&(position[axis][obj_id]), _mm256_fmadd_pd(_mm256_load_pd(&(velocity[axis][obj_id])), cps, _mm256_load_pd(&(position[axis][obj_id]))));
					}
				}
			}
		}
		barrier->arrive_and_wait();
		std::call_once(sync_flags[sync_flag_index++], [&] {
			tc.stop();
			*event_log_buffer << period_ord << "	" << tc.measured_timespan().count() << _endl_;
		});
		barrier->arrive_and_wait();
	}
	{ //after-sim log
		std::call_once(sync_flags[sync_flag_index++], [&] {
			event_log_buffer = std::make_unique<PCL::C_Event_Log_Buffer>(event_log_obj, false, "THREAD LOG");
		});
		barrier->arrive_and_wait();
		{
			std::unique_lock<std::mutex> lck(mut_ex);
			*event_log_buffer << "Thread " << thr_index << " finished working!" << _endl_;
		}
		barrier->arrive_and_wait();
		std::call_once(sync_flags[sync_flag_index++], [&] {
			event_log_buffer.reset();
			time_counter.stop();
		});
		if (settings_obj->get_settings().calculate_energy_error) {
			{
				auto temp_energy = system_energy();
				std::unique_lock<std::mutex> lck(mut_ex);
				E_c_k += temp_energy;
			}
			barrier->arrive_and_wait();
			std::call_once(sync_flags[sync_flag_index++], [&] {
				event_log_obj() << nm(E_c_k) << _endl_;
			});
		}
	} //after-sim log end
}

inline void C_Universe::data_log() {
	PCL::C_Event_Log_Buffer{ event_log_obj, false, "SIM_BASIC_DATA" }
		<< "sim_basic_data().num_of_objects = " << sim_basic_data().num_of_objects << _endl_
		<< "sim_basic_data().sim_duration = " << sim_basic_data().sim_duration << _endl_
		<< "sim_basic_data().cycles_per_second_exponent = " << sim_basic_data().cycles_per_second_exponent << _endl_
		<< "period_lenght = " << period_lenght << _endl_
		<< "periods = " << periods << _endl_;
}

inline void C_Universe::data_display() {
	std::cout << fixed
		<< line << '\n'
		<< "Ilosc obiektow: " << sim_basic_data().num_of_objects << '\n'
		<< "Czas trwania symulacji: " << sim_basic_data().sim_duration << '\n'
		<< "Ilosc cykli na sekunde: " << sim_basic_data().cycles_per_second_exponent << '\n';
	for (auto& obj : obj_arr) {
		//cout << obj;
	}
	std::cout << scientific;
}

inline void C_Universe::data_main_input() {
	do {
		if (sim_basic_data().num_of_objects < 2) {
			std::cout << line << '\n'
				<< "Minimalna ilosc obiektow to 2!" << '\n';
		}
		std::cout << line << '\n'
			<< "Podaj liczbe obiektow: ";
		cin >> sim_basic_data().num_of_objects;
	} while (sim_basic_data().num_of_objects < 2 || input_error());
	do {
		if (sim_basic_data().sim_duration == 0) {
			std::cout << line << '\n'
				<< "Czas trwania symulacji nie moze byc zerowy!" << '\n';
		}
		std::cout << line << '\n'
			<< "Podaj czas trwania symulacji (w okresach po 131072 sekund): ";
		cin >> sim_basic_data().sim_duration;
	} while (sim_basic_data().sim_duration == 0 || input_error());
	do {
		std::cout << line << '\n'
			<< "Podaj liczbe cykli na sekunde: ";
		cin >> sim_basic_data().cycles_per_second_exponent;
	} while (input_error());
}

inline void C_Universe::data_change() {
	do {
		do {
			std::cout << line << '\n'
				<< "0 - Kontynuuj bez zmian" << '\n'
				<< "1 - Czas trwania symulacji: " << sim_basic_data().sim_duration << '\n'
				<< "2 - Wykladnik cykli cykli na sekunde: " << sim_basic_data().cycles_per_second_exponent << '\n'
				<< "3 - Zmiana danych obiektu" << '\n'
				<< line << '\n'
				<< "Wybierz parametr do zmiany: ";
			cin >> mod_choice;
		} while (input_error() || incorrect_value(mod_choice > 3));
		event_log_obj() << "mod_choice = " << mod_choice << _endl_;
		switch (mod_choice) {
		case 0: {
			break;
		}
		case 1:
		{
			do {
				std::cout << line << '\n'
					<< "Podaj czas trwania symulacji w sekundach: ";
				cin >> sim_basic_data().sim_duration;
			} while (incorrect_value(sim_basic_data().sim_duration == 0) || input_error());
			break;
		}
		case 2:
		{
			do {
				std::cout << line << '\n'
					<< "Podaj wykladnik liczby cykli na sekunde: ";
				cin >> sim_basic_data().cycles_per_second_exponent;
			} while (input_error());
			break;
		}
		case 3: {
			uint_fast64_t object_number = 0;
			do {
				std::cout << line << '\n'
					<< "Numer obiektu do edytowania: ";
				cin >> object_number;
			} while (input_error() || incorrect_value(object_number >= sim_basic_data().num_of_objects));
			break;
		}
		}
	} while (mod_choice != 0);
}

inline void C_Universe::data_speed_test() {
	sim_basic_data() = settings_obj->get_settings().test_sim_basic_data();
}

inline void C_Universe::result_recording() {
	if (include_analysis) {
		bool files_good = true;
		PCL::C_Event_Log_Buffer buffer(event_log_obj, true, "OUTPUT_FILES");
		filesystem::path file_path = "simulations/simulation[" + to_string(settings_obj->get_settings().num_of_sim) + "]" + string(action_identifier == 4 ? "[ST]" : "") + "/result.bin";
		fstream file(file_path, ios::out | ios::binary);
		if (!(files_good &= file.is_open() && file)) {
			buffer() << "Creating file " << file_path.string() << "- failure" << _endl_;
		}
		buffer() << string("Result recording ") << operation_evaluation<false>(files_good) << _endl_;
	}
}

inline void C_Universe::simulation_preperator() {
	for_each(execution::par_unseq, obj_arr.begin(), obj_arr.end(), [](C_Object<dimensions>& obj)->void {
		obj.calculate_sgp();
	});
	const auto n__ = 20;
	const auto num_of_available_threads = (settings_obj->get_settings().num_of_threads_override == 0 ? thread::hardware_concurrency() : settings_obj->get_settings().num_of_threads_override);
	num_of_threads_in_use = clamp<uint_fast64_t>(floor(double(1) / n__ * sim_basic_data().num_of_objects + (n__ - 1) / (2 * n__)), 1, num_of_available_threads);
	event_log_obj() << "Number of threads - " << num_of_threads_in_use << _endl_;
	barrier = std::make_unique<PCL::C_Barrier>(num_of_threads_in_use);
	indexer = std::make_unique<PCL::C_Indexer<>>(num_of_threads_in_use);
	auto thr_binded = bind(&C_Universe::simulation_packed_thread, this);
	PCL::C_thread_set().run(thr_binded, num_of_threads_in_use);
	event_log_obj() << "Simulation's threads synchronized!" << _endl_;
}