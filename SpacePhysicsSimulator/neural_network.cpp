#include "neural_network.h"
#include <iostream>

std::vector<double> conv_vec(Eigen::VectorXd vec) {
	return { vec.data(), vec.data() + vec.rows() };
}

Eigen::VectorXd conv_vec(const std::vector<double> vec) {
	return Eigen::Map<const Eigen::VectorXd, Eigen::Unaligned>(vec.data(), vec.size());
}

void neural_network::calc_next_layer(Eigen::VectorXd& data, size_t layer_index) const {
	data = (weights[layer_index] * data - biases[layer_index]).unaryExpr(sigmoid);
}

std::vector<double> neural_network::calc_result(const std::vector<double> data) const {
	Eigen::VectorXd result = conv_vec(data);
	for (size_t i = 0; i < layer_sizes.size() - 1; i++) {
		calc_next_layer(result, i);
	}
	return conv_vec(result);
}

std::pair<Eigen::VectorXd, Eigen::VectorXd> neural_network::calc_next_layer_with_z(Eigen::VectorXd data, size_t layer_index) {
	auto z_vec = weights[layer_index] * data - biases[layer_index];
	auto a_vec = z_vec.unaryExpr(sigmoid);
	return { a_vec, z_vec };
}

std::pair<std::vector<Eigen::VectorXd>, std::vector<Eigen::VectorXd>> neural_network::calc_result_with_z(Eigen::VectorXd data) {
	std::pair<std::vector<Eigen::VectorXd>, std::vector<Eigen::VectorXd>> results = { {data},{{}} };
	for (size_t i = 0; i < layer_sizes.size() - 1; i++) {
		auto [a, z] = calc_next_layer_with_z(results.first.back(), i);
		results.first.push_back(a);
		results.second.push_back(z);
	}
	return results;
}

Eigen::VectorXd neural_network::calc_result_internal(Eigen::VectorXd data) const {
	for (size_t i = 0; i < layer_sizes.size() - 1; i++) {
		calc_next_layer(data, i);
	}
	return data;
}

neural_network::neural_network(std::vector<size_t> layer_sizes_arg) : layer_sizes(layer_sizes_arg) {
	weights.resize(layer_sizes.size() - 1);
	weights_err.resize(layer_sizes.size() - 1);
	biases.resize(layer_sizes.size() - 1);
	biases_err.resize(layer_sizes.size() - 1);
	for (size_t i = 0; i < weights.size(); i++) {
		weights[i] = Eigen::MatrixXd::Random(layer_sizes[i + 1], layer_sizes[i]);
		weights_err[i] = Eigen::MatrixXd::Zero(layer_sizes[i + 1], layer_sizes[i]);
		biases[i] = Eigen::VectorXd::Zero(layer_sizes[i + 1]);
		biases_err[i] = Eigen::VectorXd::Zero(layer_sizes[i + 1]);
	}
}

void neural_network::train(const std::vector<double>& data, std::vector<double> expected_result) {
	auto er = conv_vec(expected_result);
	auto data_e = conv_vec(data);
	auto [a, z] = calc_result_with_z(data_e);
	auto diff = 2 * (er - a.back());
	biases_err.back() = diff * z.back().unaryExpr(sigmoid_derivate);
	for (size_t i = layer_sizes.size() - 2; i > 0; i--) { //calculate delta/bias error
		biases_err[i - 1] = (weights[i].transpose() * biases_err[i]) * z[i].unaryExpr(sigmoid_derivate);
	}
	for (size_t i = layer_sizes.size() - 2; i <= layer_sizes.size() - 2; i--) { //calculate weight error
		weights_err[i] = biases_err[i] * a[i].transpose();
	}
	std::transform(std::execution::par_unseq, weights.begin(), weights.end(), weights_err.begin(), weights.begin(), [&](Eigen::MatrixXd& weights_one_layer, Eigen::MatrixXd& weights_err_one_layer) {
		return weights_one_layer + weights_err_one_layer * learning_coef;
	});
	std::transform(std::execution::par_unseq, biases.begin(), biases.end(), biases_err.begin(), biases.begin(), [&](Eigen::VectorXd& biases_one_layer, Eigen::VectorXd& biases_err_one_layer) {
		return biases_one_layer + biases_err_one_layer * learning_coef;
	});
}

void neural_network::print_final() {
	for (auto& mt : weights) {
		std::cout << mt << '\n';
	}
	for (auto& b : biases) {
		std::cout << b << '\n';
	}
}
