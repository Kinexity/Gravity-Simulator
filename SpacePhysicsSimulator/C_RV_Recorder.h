#pragma once
#ifndef C_RV_R
#define C_RV_R
#include <tuple>
#include <map>
#include <functional>

template <typename Ret_Type, typename... Func_Args>
class C_RV_Storage {
private:
	std::map<std::tuple<Func_Args...>, Ret_Type> ret_val_map;
	std::function<Ret_Type(Func_Args...)> func_ptr;
public:
	C_RV_Storage(std::function<Ret_Type(Func_Args...)>& f_ptr) : func_ptr(f_ptr) {};
	C_RV_Storage() = default;
	~C_RV_Storage() = default;
	decltype(func_ptr)& operator=(std::function < Ret_Type(Func_Args...)> f_ptr) {
		return (func_ptr = f_ptr);
	}
	Ret_Type operator()(Func_Args... args) {
		const auto it = ret_val_map.find({ args... });
		return (it != ret_val_map.end() ? it->second : (ret_val_map[{ args... }] = func_ptr(args...)));
	};
};

#endif // !C_RV_R