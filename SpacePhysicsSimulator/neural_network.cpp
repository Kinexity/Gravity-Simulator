#include "neural_network.h"
#include <iostream>

std::vector<double> conv_vec(Eigen::VectorXd vec) {
	return { vec.data(), vec.data() + vec.rows() };
}

Eigen::VectorXd conv_vec(const std::vector<double> vec) {
	return Eigen::Map<const Eigen::VectorXd, Eigen::Unaligned>(vec.data(), vec.size());
}

void neural_network::calc_next_layer(Eigen::VectorXd& data, size_t layer_index) {
	Eigen::VectorXd next_data;
	auto& current_weights = weights[layer_index];
	auto& current_biases = biases[layer_index];
	next_data = (current_weights * data - current_biases).unaryExpr(sigmoid);
	data = next_data;
}

std::vector<double> neural_network::calc_result(const std::vector<double> data) {
	Eigen::VectorXd result = conv_vec(data);
	for (size_t i = 0; i < layer_sizes.size() - 1; i++) {
		calc_next_layer(result, i);
	}
	return conv_vec(result);
}

std::pair<Eigen::VectorXd, Eigen::VectorXd> neural_network::calc_next_layer_with_z(Eigen::VectorXd data, size_t layer_index) {
	auto& current_weights = weights[layer_index];
	auto& current_biases = biases[layer_index];
	auto z_vec = current_weights * data - current_biases;
	auto a_vec = z_vec.unaryExpr(sigmoid);
	return { a_vec, z_vec };
}

std::vector<std::pair<Eigen::VectorXd, Eigen::VectorXd>> neural_network::calc_result_with_z(const std::vector<double> data) {
	std::vector<std::pair<Eigen::VectorXd, Eigen::VectorXd>> results;
	results.push_back({ conv_vec(data), {} });
	for (size_t i = 0; i < layer_sizes.size() - 1; i++) {
		results.push_back(calc_next_layer_with_z(results.back().first, i));
	}
	return results;
}

neural_network::neural_network(std::vector<size_t> layer_sizes_arg) : layer_sizes(layer_sizes_arg) {
	dummy_vector.resize(*std::max_element(layer_sizes.begin(), layer_sizes.end()));
	weights.resize(layer_sizes.size() - 1);
	for (size_t i = 0; i < weights.size(); i++) {
		weights[i] = Eigen::MatrixXd::Random(layer_sizes[i + 1], layer_sizes[i]);
	}
	weights_err = weights;
	biases.resize(layer_sizes.size() - 1);
	for (size_t i = 0; i < weights.size(); i++) {
		biases[i] = Eigen::VectorXd::Random(layer_sizes[i + 1]);
	}
	biases_err = biases;
}

void neural_network::train(const std::vector<double>& data, std::vector<double> expected_result) {
	auto er = conv_vec(expected_result);
	auto results = calc_result_with_z(data);
	Eigen::VectorXd
		calculated_result = results.back().first,
		diff,
		intermediate_res;
	diff = 2 * (er - calculated_result);
	biases_err.back() = diff * results.back().second.unaryExpr(sigmoid_derivate);
	for (size_t i = layer_sizes.size() - 2; i > 0; i--) { //calculate delta/bias error
		Eigen::MatrixXd weights_layer_copy = weights[i];
		Eigen::VectorXd& biases_err_local = biases_err[i - 1];
		Eigen::VectorXd& z_local = results[i].second;
		biases_err_local = (Eigen::Transpose(weights_layer_copy)* biases_err[i])* z_local.unaryExpr(sigmoid_derivate);
	}
	for (size_t i = layer_sizes.size() - 2; i <= layer_sizes.size() - 2; i--) { //calculate weight error
		auto& weights_err_layer = weights_err[i];
		auto& deltas = biases_err[i];
		auto& a = results[i].first;
		weights_err_layer = deltas * Eigen::Transpose(a);
	}
	std::transform(std::execution::par_unseq, weights.begin(), weights.end(), weights_err.begin(), dummy_vector.begin(), [&](Eigen::MatrixXd& weights_one_layer, Eigen::MatrixXd& weights_err_one_layer) {
		weights_one_layer += weights_err_one_layer * learning_coef;
		return 0;
	});
	std::transform(std::execution::par_unseq, biases.begin(), biases.end(), biases_err.begin(), dummy_vector.begin(), [&](Eigen::VectorXd& biases_one_layer, Eigen::VectorXd& biases_err_one_layer) {
		biases_one_layer += biases_err_one_layer * learning_coef;
		return 0;
	});
}