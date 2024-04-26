#pragma once
#ifndef C_RV_R
#define C_RV_R
#include <tuple>
#include <map>
#include <set>
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
		clear();
		return (func_ptr = f_ptr);
	}
	Ret_Type operator()(Func_Args... args) {
		auto packed = std::make_tuple(args...);
		const auto it = ret_val_map.find(packed);
		return (it != ret_val_map.end() ? it->second : (ret_val_map[packed] = func_ptr(args...)));
	};
	void clear() {
		ret_val_map.clear();
	}
};

template <typename Ret_Type, typename... Func_Args>
class return_value_mem_eff_wrapper {
private:
	std::map<std::tuple<Func_Args...>, typename std::set<Ret_Type>::iterator> ret_val_map;
	std::set<Ret_Type> ret_val_set;
	std::function<Ret_Type(Func_Args...)> func_ptr;
public:
	return_value_mem_eff_wrapper(std::function<Ret_Type(Func_Args...)>& f_ptr) : func_ptr(f_ptr) {};
	return_value_mem_eff_wrapper() = default;
	~return_value_mem_eff_wrapper() = default;
	decltype(func_ptr)& operator=(std::function < Ret_Type(Func_Args...)> f_ptr) {
		clear();
		return (func_ptr = f_ptr);
	}
	Ret_Type operator()(Func_Args... args) {
		auto packed = std::make_tuple(args...);
		const auto it = ret_val_map.find(packed);
		if (it != ret_val_map.end()) {
			return *(it->second);
		}
		else {
			auto [set_it, success] = ret_val_set.insert(func_ptr(args...));
			ret_val_map[packed] = set_it;
			return *set_it;
		}
	};
	void clear() {
		ret_val_map.clear();
		ret_val_set.clear();
	}
};

template <typename Ret_Type, typename... Func_Args>
class atomic_return_value_wrapper {
private:
	std::map<std::tuple<typename std::remove_cvref<Func_Args>::type...>, Ret_Type> ret_val_map;
	std::function<Ret_Type(Func_Args...)> func_ptr;
	std::mutex sync_mtx;
public:
	atomic_return_value_wrapper(std::function<Ret_Type(Func_Args...)>& f_ptr) : func_ptr(f_ptr) {};
	atomic_return_value_wrapper() = default;
	~atomic_return_value_wrapper() = default;
	decltype(func_ptr)& operator=(std::function < Ret_Type(Func_Args...)> f_ptr) {
		std::lock_guard<std::mutex> lock(sync_mtx);
		clear();
		return (func_ptr = f_ptr);
	}
	Ret_Type operator()(Func_Args... args) {
		auto packed = std::make_tuple(args...);
		typename decltype(ret_val_map)::iterator it;
		{
			std::lock_guard<std::mutex> lock(sync_mtx);
			it = ret_val_map.find(packed);
		}
		Ret_Type return_value;
		if (it == ret_val_map.end()) {
			return_value = func_ptr(args...);
			std::lock_guard<std::mutex> lock(sync_mtx);
			ret_val_map[packed] = return_value;
		}
		else {
			return_value = it->second;
		}
		return return_value;
	};
	void clear() {
		std::lock_guard<std::mutex> lock(sync_mtx);
		ret_val_map.clear();
	}
};

#endif // !C_RV_R