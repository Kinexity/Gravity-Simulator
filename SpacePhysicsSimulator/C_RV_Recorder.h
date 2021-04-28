#pragma once
#ifndef C_RV_R
#define C_RV_R
#include <tuple>
#include <map>
#include <functional>

template <typename Ret_Type, typename... Func_Args>
class return_value_wrapper {
private:
	std::map<std::tuple<Func_Args...>, Ret_Type> ret_val_map;
	std::function<Ret_Type(Func_Args...)> func_ptr;
public:
	return_value_wrapper(std::function<Ret_Type(Func_Args...)>& f_ptr) : func_ptr(f_ptr) {};
	return_value_wrapper() = default;
	~return_value_wrapper() = default;
	decltype(func_ptr)& operator=(std::function < Ret_Type(Func_Args...)> f_ptr) {
		ret_val_map.clear();
		return (func_ptr = f_ptr);
	}
	Ret_Type operator()(Func_Args... args) {
		auto packed = std::make_tuple(args...);
		const auto it = ret_val_map.find(packed);
		return (it != ret_val_map.end() ? it->second : (ret_val_map[packed] = func_ptr(args...)));
	};
};

#endif // !C_RV_R