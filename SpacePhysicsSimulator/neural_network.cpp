#include "neural_network.h"

std::vector<double>&& neural_network::calc_next_layer(std::vector<double> data, size_t layer_index) {
	std::vector<double> next_data(layer_sizes[layer_index + 1]);
	auto& current_weights = weights[layer_index];
	auto& current_biases = biases[layer_index];
	std::transform(std::execution::par_unseq, current_weights.begin(), current_weights.end(), current_biases.begin(), next_data.begin(), [&](std::vector<double>& weights_row, double& bias) {
		return sigmoid(std::transform_reduce(std::execution::par_unseq, data.begin(), data.end(), weights_row.begin(), 0., std::plus<>(), std::multiplies<>()) - bias);
	});
	return std::move(next_data);
}

std::vector<double>&& neural_network::calc_result(const std::vector<double> data) {
	std::vector<double> result = data;
	for (size_t i = 0; i < layer_sizes.size() - 1; i++) {
		result = calc_next_layer(result, i);
	}
	return std::move(result);
}

neural_network::neural_network(std::vector<size_t> layer_sizes_arg) : layer_sizes(layer_sizes_arg) {
	dummy_vector.resize(*std::max_element(layer_sizes.begin(), layer_sizes.end()));
	weights.resize(layer_sizes.size() - 1);
	std::for_each(std::execution::par_unseq, weights.begin(), weights.end(), [&](std::vector<std::vector<double>>& weights_one_layer) {
		auto one_layer_index = std::distance(&*weights.begin(), &weights_one_layer);
		weights_one_layer.resize(layer_sizes[one_layer_index + 1]);
		std::for_each(std::execution::par_unseq, weights_one_layer.begin(), weights_one_layer.end(), [&](std::vector<double>& weights_one_row) {
			weights_one_row.resize(layer_sizes[one_layer_index]);
			std::generate(std::execution::par_unseq, weights_one_row.begin(), weights_one_row.end(), [&]() {return uni_real_rand(rd); });
		});
	});
	weights_err = weights;
	biases.resize(layer_sizes.size() - 1);
	std::for_each(std::execution::par_unseq, biases.begin(), biases.end(), [&](std::vector<double>& biases_one_layer) {
		auto one_layer_index = std::distance(&*biases.begin(), &biases_one_layer);
		biases_one_layer.resize(layer_sizes[one_layer_index + 1]);
		std::generate(std::execution::par_unseq, biases_one_layer.begin(), biases_one_layer.end(), [&]() {return uni_real_rand(rd); });
	});
	biases_err = biases;
}

void neural_network::train(const std::vector<double>& data, std::vector<double> expected_result) {

	std::transform(std::execution::par_unseq, weights.begin(), weights.end(), weights_err.begin(), dummy_vector.begin(), [&](std::vector<std::vector<double>>& weights_one_layer, std::vector<std::vector<double>>& weights_err_one_layer) {
		std::transform(std::execution::par_unseq, weights_one_layer.begin(), weights_one_layer.end(), weights_err_one_layer.begin(), dummy_vector.begin(), [&](std::vector<double>& weights_one_row, std::vector<double>& weights_err_one_row) {
			std::transform(std::execution::par_unseq, weights_one_row.begin(), weights_one_row.end(), weights_err_one_row.begin(), weights_one_row.begin(), [&](double& weight, double& weight_err) {return std::fma(weight_err, learning_coef, weight); });
			return 0;
		});
		return 0;
	});
	std::transform(std::execution::par_unseq, biases.begin(), biases.end(), biases_err.begin(), dummy_vector.begin(), [&](std::vector<double>& biases_one_layer, std::vector<double>& biases_err_one_layer) {
		std::transform(std::execution::par_unseq, biases_one_layer.begin(), biases_one_layer.end(), biases_err_one_layer.begin(), biases_one_layer.begin(), [&](double& bias, double& bias_err) {return std::fma(bias_err, learning_coef, bias); });
		return 0;
	});
}