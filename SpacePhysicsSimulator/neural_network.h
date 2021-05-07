#pragma once
#ifndef neural_network_h
#define neural_network_h
#include <vector>
#include <intrin.h>
#include <algorithm>
#include <numeric>
#include <execution>
#include <functional>
#include <array>
#include <numeric>
#include "PureCPPLib/C_Random.h"

class neural_network {
private:
	std::vector<size_t>
		layer_sizes,
		dummy_vector; //1 + deep + 1
	std::vector<std::vector<std::vector<double>>>
		weights,
		weights_err;
	std::vector<std::vector<double>>
		biases,
		biases_err; //is equal to neuron error
	std::random_device
		rd;
	double
		learning_coef = 0.01;
	std::uniform_real_distribution<double>
		uni_real_rand{ -10.,10. };
	std::function<double(double)>
		sigmoid = [&](double x) { return (std::erf(x) + 1) / 2; },
		sigmoid_derivate = [&](double x) { return std::exp(-x * x) / 1.7724538509055160273; };
	void
		calc_next_layer(std::vector<double>& data, size_t layer_index);
	std::pair<std::vector<double>, std::vector<double>>
		calc_next_layer_with_z(std::vector<double> data, size_t layer_index);
	std::vector<std::pair<std::vector<double>, std::vector<double>>>
		calc_result_with_z(std::vector<double> data);
public:
	neural_network(std::vector<size_t> layer_sizes_arg);
	~neural_network() = default;
	void
		train(const std::vector<double>& data, std::vector<double> expected_result);
	std::vector<double>
		calc_result(std::vector<double> data);
};

#endif // !neural_network_h