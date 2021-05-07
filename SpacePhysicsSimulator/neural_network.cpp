#include "neural_network.h"
#include <iostream>

void neural_network::calc_next_layer(std::vector<double>& data, size_t layer_index) {
	std::vector<double> next_data(layer_sizes[layer_index + 1]);
	auto& current_weights = weights[layer_index];
	auto& current_biases = biases[layer_index];
	std::transform(std::execution::par_unseq, current_weights.begin(), current_weights.end(), current_biases.begin(), next_data.begin(), [&](std::vector<double>& weights_row, double& bias) {
		return sigmoid(std::transform_reduce(std::execution::par_unseq, data.begin(), data.end(), weights_row.begin(), 0., std::plus<>(), std::multiplies<>()) - bias);
	});
	//std::cout << data.size() << "\n";
	data = next_data;
}

std::vector<double> neural_network::calc_result(const std::vector<double> data) {
	std::vector<double> result = data;
	for (size_t i = 0; i < layer_sizes.size() - 1; i++) {
		calc_next_layer(result, i);
		//std::cout << result.size() << "\n";
	}
	return result;
}

std::pair<std::vector<double>, std::vector<double>> neural_network::calc_next_layer_with_z(std::vector<double> data, size_t layer_index) {
	std::vector<double> next_data(layer_sizes[layer_index + 1]);
	auto& current_weights = weights[layer_index];
	auto& current_biases = biases[layer_index];
	std::transform(std::execution::par_unseq, current_weights.begin(), current_weights.end(), current_biases.begin(), next_data.begin(), [&](std::vector<double>& weights_row, double& bias) {
		return std::transform_reduce(std::execution::par_unseq, data.begin(), data.end(), weights_row.begin(), 0., std::plus<>(), std::multiplies<>()) - bias;
	});
	std::vector<double> z_vec = next_data;
	std::transform(std::execution::par_unseq, next_data.begin(), next_data.end(), next_data.begin(), sigmoid);
	return std::move(std::pair<std::vector<double>, std::vector<double>>{ next_data, z_vec});
}

std::vector<std::pair<std::vector<double>, std::vector<double>>> neural_network::calc_result_with_z(const std::vector<double> data) {
	std::vector<std::pair<std::vector<double>, std::vector<double>>> results;
	results.push_back({ data, std::vector<double>() });
	for (size_t i = 0; i < layer_sizes.size() - 1; i++) {
		results.push_back(calc_next_layer_with_z(results.back().first, i));
	}
	return std::move(results);
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

std::vector<double> operator+(std::vector<double> v1, std::vector<double> v2) {
	std::vector<double> res(v1.size());
	std::transform(std::execution::par_unseq, v1.begin(), v1.end(), v2.begin(), res.begin(), std::plus());
	return res;
}

void neural_network::train(const std::vector<double>& data, std::vector<double> expected_result) {
	auto results = calc_result_with_z(data);
	std::vector<double>
		calculated_result = results.back().first,
		diff(results.back().first.size()),
		intermediate_res;
	std::transform(std::execution::par_unseq, calculated_result.begin(), calculated_result.end(), expected_result.begin(), diff.begin(), [&](double cr, double er) { return 2 * (er - cr); });
	std::transform(std::execution::par_unseq, diff.begin(), diff.end(), results.back().second.begin(), biases_err.back().begin(), [&](double diff_C, double z_L) { return diff_C * sigmoid_derivate(z_L); });
	for (size_t i = layer_sizes.size() - 2; i > 0; i--) { //calculate delta/bias error
		std::vector<std::vector<double>> weights_layer_copy = weights[i];
		std::vector<double>& biases_err_local = biases_err[i - 1];
		std::vector<double>& z_local = results[i].second;
		biases_err_local = std::transform_reduce(std::execution::par_unseq, weights_layer_copy.begin(), weights_layer_copy.end(), biases_err[i].begin(), std::vector<double>(biases_err[i - 1].size()), [&](std::vector<double> v1, std::vector<double> v2) { return v1 + v2; }, [&](std::vector<double>& row, double bias) {
			std::for_each(std::execution::par_unseq, row.begin(), row.end(), [&](double& elem) {
				elem *= bias;
			});
			return row;
		});
		std::transform(std::execution::par_unseq, biases_err_local.begin(), biases_err_local.end(), z_local.begin(), biases_err_local.begin(), [&](double bias, double z) { return bias * sigmoid_derivate(z); });
	}
	for (size_t i = layer_sizes.size() - 2; i <= layer_sizes.size() - 2; i--) { //calculate weight error
		std::vector<std::vector<double>>& weights_err_layer = weights_err[i];
		std::vector<double>& deltas = biases_err[i];
		std::vector<double>& a = results[i].first;
		std::transform(std::execution::par_unseq, weights_err_layer.begin(), weights_err_layer.end(), deltas.begin(), dummy_vector.begin(), [&](std::vector<double>& row, double delta) {
			std::transform(std::execution::par_unseq, a.begin(), a.end(), row.begin(), [&](double elem) { return elem * delta; });
			return 0;
		});
	}
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