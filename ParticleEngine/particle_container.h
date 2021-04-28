#pragma once
#ifndef particle_container_h
#define particle_container_h
#include <vector>
#include <array>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <execution>
#include "PureCPPLib/C_Random.h"

template<size_t dimensions>
class particle_container {
private:
	PCL::C_Random
		rand;
	std::array<std::vector<double>, dimensions>
		position,
		velocity;
	double
		particle_radius,
		particle_mass,
		container_side_half_lenght,
		initial_energy;
	size_t
		number_of_sectors,
		sectors_per_side,
		number_of_particles;
	double
		total_energy(),
		particle_distance(size_t particle_index_1, size_t particle_index_2),
		speed_difference_distance(size_t particle_index_1, size_t particle_index_2);
	std::vector<std::vector<size_t>>
		particles_in_sector;
public:
	particle_container() = default;
	~particle_container() = default;
	void
		simulation(),
		simulate_sector(size_t sector_index),
		randomise_positions();
};

template<size_t dimensions>
inline double particle_container<dimensions>::total_energy() {
	return particle_mass * std::transform_reduce(std::execution::par_unseq, velocity.begin(), velocity.end(), 0., std::plus<>(), [&](std::vector<double>& v_) {
		return std::transform_reduce(std::execution::par_unseq, v_.begin(), v_.end(), 0., std::plus<>(), [&](double& vel) { return vel * vel; });
	});
}

template<size_t dimensions>
inline double particle_container<dimensions>::speed_difference_distance(size_t particle_index_1, size_t particle_index_2)
{
	return 0.0;
}

template<size_t dimensions>
inline double particle_container<dimensions>::particle_distance(size_t particle_index_1, size_t particle_index_2) {
	auto distance_2 = 0.;
	for (size_t dimension_index = 0; dimension_index < dimensions; dimension_index++) {
		auto x_diff = position[dimension_index][particle_index_1] - position[dimension_index][particle_index_2];
		distance_2 += x_diff * x_diff;
	}
	return std::sqrt(distance_2);
}

template<size_t dimensions>
inline void particle_container<dimensions>::simulation() {
	simulate_sector(0);
}

template<size_t dimensions>
inline void particle_container<dimensions>::randomise_positions() {
}

template<size_t dimensions>
inline void particle_container<dimensions>::simulate_sector(size_t sector_index) {
	for (auto& particle_index : particles_in_sector[sector_index]) {

	}
	for (size_t sector_index = 0; sector_index < sector_index + number_of_sectors; sector_index *= number_of_sectors) {
		simulate_sector(sector_index + number_of_sectors);
	}
}

#endif particle_container_h