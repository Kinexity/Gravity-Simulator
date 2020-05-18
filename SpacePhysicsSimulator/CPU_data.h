#pragma once
#include <iostream>  
#include <vector>  
#include <bitset>  
#include <array>  
#include <string>  
#include <intrin.h> 
using namespace std;

class InstructionSet_internal {
public:
	InstructionSet_internal();
	int_fast32_t
		nIds_,
		nExIds_;
	string
		vendor_,
		brand_;
	bool
		isintel_,
		isAMD_;
	bitset<32>
		f_1_ECX_,
		f_1_EDX_,
		f_7_EBX_,
		f_7_ECX_,
		f_81_ECX_,
		f_81_EDX_;
	vector<array<int_fast32_t, 4>>
		data_,
		extdata_;
};

class InstructionSet {
public:
	static string
		Vendor(),
		Brand();
	static bool
		SSE3(),
		PCLMULQDQ(),
		MONITOR(),
		SSSE3(),
		FMA(),
		CMPXCHG16B(),
		SSE41(),
		SSE42(),
		MOVBE(),
		POPCNT(),
		AES(),
		XSAVE(),
		OSXSAVE(),
		AVX(),
		F16C(),
		RDRAND(),
		MSR(),
		CX8(),
		SEP(),
		CMOV(),
		CLFSH(),
		MMX(),
		FXSR(),
		SSE(),
		SSE2(),
		FSGSBASE(),
		BMI1(),
		HLE(),
		AVX2(),
		BMI2(),
		ERMS(),
		INVPCID(),
		RTM(),
		AVX512F(),
		RDSEED(),
		ADX(),
		AVX512PF(),
		AVX512ER(),
		AVX512CD(),
		SHA(),
		PREFETCHWT1(),
		LAHF(),
		LZCNT(),
		ABM(),
		SSE4a(),
		XOP(),
		TBM(),
		SYSCALL(),
		MMXEXT(),
		RDTSCP(),
		_3DNOWEXT(),
		_3DNOW();
	static InstructionSet_internal CPU_Rep;
};

string main_CPU();

void additional_CPU_info();