#pragma once
#ifndef C_Tt
#define C_Tt
#define _ENABLE_EXTENDED_ALIGNED_STORAGE
#include <string>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <algorithm>
#include <memory>
#include <random>
#include <chrono>
#include <list>
#include <tuple>
#include <bitset>
#include <vector>
#include <iterator>
#include <array>
#include <condition_variable>
#include <conio.h>
#include <atomic>
#include <intrin.h>
#include <map>
#include <functional>
#include <ranges>
#include "C_RV_Recorder.h"
#include "PureCPPLib/C_Indexer.h"
#include "PureCPPLib/C_Time_Counter.h"
#include "PureCPPLib/C_Random.h"
#include "PureCPPLib/C_XML.h"
#include "constants.h"

template < std::size_t I = 0, typename...Tp >
inline void sumTuples(const std::tuple < Tp...>& t1, const std::tuple < Tp...>& t2, std::tuple < Tp...>& _result) noexcept {
	if constexpr (I < sizeof...(Tp)) {
		std::get < I >(_result) = std::get < I >(t1) + std::get < I >(t2);
		sumTuples < I + 1, Tp...>(t1, t2, _result);
	}
}

template < typename...Tp >
inline std::tuple < Tp...> operator+(const std::tuple < Tp...>& t1, const std::tuple < Tp...>& t2) {
	std::tuple < Tp...> _result;
	sumTuples < 0, Tp...>(t1, t2, _result);
	return _result;
}

template < typename...Tp >
inline std::tuple < Tp...> operator+=(std::tuple < Tp...>& t1, const std::tuple < Tp...>& t2) {
	sumTuples < 0, Tp...>(t1, t2, t1);
	return t1;
}

class C_Test {
public:
	PCL::C_Event_Log_Base&
		event_log;
	C_Test(PCL::C_Event_Log_Base& event_log_ref);
	~C_Test() = default;
	void
		run();
};

#endif // !C_Tt