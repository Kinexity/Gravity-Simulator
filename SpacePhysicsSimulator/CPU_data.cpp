#include "CPU_data.h"

InstructionSet_internal InstructionSet::CPU_Rep;

std::string InstructionSet::Vendor() { return InstructionSet::CPU_Rep.vendor_; }
std::string InstructionSet::Brand() { return InstructionSet::CPU_Rep.brand_; }

bool InstructionSet::InstructionSet::SSE3() { return InstructionSet::CPU_Rep.f_1_ECX_[0]; }
bool InstructionSet::InstructionSet::PCLMULQDQ() { return InstructionSet::CPU_Rep.f_1_ECX_[1]; }
bool InstructionSet::InstructionSet::MONITOR() { return InstructionSet::CPU_Rep.f_1_ECX_[3]; }
bool InstructionSet::InstructionSet::SSSE3() { return InstructionSet::CPU_Rep.f_1_ECX_[9]; }
bool InstructionSet::InstructionSet::FMA() { return InstructionSet::CPU_Rep.f_1_ECX_[12]; }
bool InstructionSet::CMPXCHG16B() { return InstructionSet::CPU_Rep.f_1_ECX_[13]; }
bool InstructionSet::SSE41() { return InstructionSet::CPU_Rep.f_1_ECX_[19]; }
bool InstructionSet::SSE42() { return InstructionSet::CPU_Rep.f_1_ECX_[20]; }
bool InstructionSet::MOVBE() { return InstructionSet::CPU_Rep.f_1_ECX_[22]; }
bool InstructionSet::POPCNT() { return InstructionSet::CPU_Rep.f_1_ECX_[23]; }
bool InstructionSet::AES() { return InstructionSet::CPU_Rep.f_1_ECX_[25]; }
bool InstructionSet::XSAVE() { return InstructionSet::CPU_Rep.f_1_ECX_[26]; }
bool InstructionSet::OSXSAVE() { return InstructionSet::CPU_Rep.f_1_ECX_[27]; }
bool InstructionSet::AVX() { return InstructionSet::CPU_Rep.f_1_ECX_[28]; }
bool InstructionSet::F16C() { return InstructionSet::CPU_Rep.f_1_ECX_[29]; }
bool InstructionSet::RDRAND() { return InstructionSet::CPU_Rep.f_1_ECX_[30]; }

bool InstructionSet::MSR() { return InstructionSet::CPU_Rep.f_1_EDX_[5]; }
bool InstructionSet::CX8() { return InstructionSet::CPU_Rep.f_1_EDX_[8]; }
bool InstructionSet::SEP() { return InstructionSet::CPU_Rep.f_1_EDX_[11]; }
bool InstructionSet::CMOV() { return InstructionSet::CPU_Rep.f_1_EDX_[15]; }
bool InstructionSet::CLFSH() { return InstructionSet::CPU_Rep.f_1_EDX_[19]; }
bool InstructionSet::MMX() { return InstructionSet::CPU_Rep.f_1_EDX_[23]; }
bool InstructionSet::FXSR() { return InstructionSet::CPU_Rep.f_1_EDX_[24]; }
bool InstructionSet::SSE() { return InstructionSet::CPU_Rep.f_1_EDX_[25]; }
bool InstructionSet::SSE2() { return InstructionSet::CPU_Rep.f_1_EDX_[26]; }

bool InstructionSet::FSGSBASE() { return InstructionSet::CPU_Rep.f_7_EBX_[0]; }
bool InstructionSet::BMI1() { return InstructionSet::CPU_Rep.f_7_EBX_[3]; }
bool InstructionSet::HLE() { return InstructionSet::CPU_Rep.isintel_ && InstructionSet::CPU_Rep.f_7_EBX_[4]; }
bool InstructionSet::AVX2() { return InstructionSet::CPU_Rep.f_7_EBX_[5]; }
bool InstructionSet::BMI2() { return InstructionSet::CPU_Rep.f_7_EBX_[8]; }
bool InstructionSet::ERMS() { return InstructionSet::CPU_Rep.f_7_EBX_[9]; }
bool InstructionSet::INVPCID() { return InstructionSet::CPU_Rep.f_7_EBX_[10]; }
bool InstructionSet::RTM() { return InstructionSet::CPU_Rep.isintel_ && InstructionSet::CPU_Rep.f_7_EBX_[11]; }
bool InstructionSet::AVX512F() { return InstructionSet::CPU_Rep.f_7_EBX_[16]; }
bool InstructionSet::RDSEED() { return InstructionSet::CPU_Rep.f_7_EBX_[18]; }
bool InstructionSet::ADX() { return InstructionSet::CPU_Rep.f_7_EBX_[19]; }
bool InstructionSet::AVX512PF() { return InstructionSet::CPU_Rep.f_7_EBX_[26]; }
bool InstructionSet::AVX512ER() { return InstructionSet::CPU_Rep.f_7_EBX_[27]; }
bool InstructionSet::AVX512CD() { return InstructionSet::CPU_Rep.f_7_EBX_[28]; }
bool InstructionSet::SHA() { return InstructionSet::CPU_Rep.f_7_EBX_[29]; }

bool InstructionSet::PREFETCHWT1() { return InstructionSet::CPU_Rep.f_7_ECX_[0]; }

bool InstructionSet::LAHF() { return InstructionSet::CPU_Rep.f_81_ECX_[0]; }
bool InstructionSet::LZCNT() { return InstructionSet::CPU_Rep.isintel_ && InstructionSet::CPU_Rep.f_81_ECX_[5]; }
bool InstructionSet::ABM() { return InstructionSet::CPU_Rep.isAMD_ && InstructionSet::CPU_Rep.f_81_ECX_[5]; }
bool InstructionSet::SSE4a() { return InstructionSet::CPU_Rep.isAMD_ && InstructionSet::CPU_Rep.f_81_ECX_[6]; }
bool InstructionSet::XOP() { return InstructionSet::CPU_Rep.isAMD_ && InstructionSet::CPU_Rep.f_81_ECX_[11]; }
bool InstructionSet::TBM() { return InstructionSet::CPU_Rep.isAMD_ && InstructionSet::CPU_Rep.f_81_ECX_[21]; }

bool InstructionSet::SYSCALL() { return InstructionSet::CPU_Rep.isintel_ && InstructionSet::CPU_Rep.f_81_EDX_[11]; }
bool InstructionSet::MMXEXT() { return InstructionSet::CPU_Rep.isAMD_ && InstructionSet::CPU_Rep.f_81_EDX_[22]; }
bool InstructionSet::RDTSCP() { return InstructionSet::CPU_Rep.isintel_ && InstructionSet::CPU_Rep.f_81_EDX_[27]; }
bool InstructionSet::_3DNOWEXT() { return InstructionSet::CPU_Rep.isAMD_ && InstructionSet::CPU_Rep.f_81_EDX_[30]; }
bool InstructionSet::_3DNOW() { return InstructionSet::CPU_Rep.isAMD_ && InstructionSet::CPU_Rep.f_81_EDX_[31]; }

InstructionSet_internal::InstructionSet_internal()
	: nIds_(0),
	nExIds_(0),
	isintel_(false),
	isAMD_(false),
	f_1_ECX_(0),
	f_1_EDX_(0),
	f_7_EBX_(0),
	f_7_ECX_(0),
	f_81_ECX_(0),
	f_81_EDX_(0),
	data_(),
	extdata_() {
	//int_fast32_t cpuInfo[4] = {-1};
	std::array<int_fast32_t, 4> cpui;
	// Calling __cpuid with 0x0 as the function_id argument
	// gets the number of the highest valid function ID.
	__cpuid(cpui.data(), 0);
	nIds_ = cpui[0];
	for (int_fast32_t i = 0; i <= nIds_; ++i) {
		__cpuidex(cpui.data(), i, 0);
		data_.push_back(cpui);
	}
	// Capture vendor std::string
	char vendor[0x20];
	memset(vendor, 0, sizeof(vendor));
	*reinterpret_cast<int_fast32_t*>(vendor) = data_[0][1];
	*reinterpret_cast<int_fast32_t*>(vendor + 4) = data_[0][3];
	*reinterpret_cast<int_fast32_t*>(vendor + 8) = data_[0][2];
	vendor_ = vendor;
	if (vendor_ == "Genuineintel") {
		isintel_ = true;
	}
	else if (vendor_ == "AuthenticAMD") {
		isAMD_ = true;
	}
	// load bitset with flags for function 0x00000001
	if (nIds_ >= 1) {
		f_1_ECX_ = data_[1][2];
		f_1_EDX_ = data_[1][3];
	}
	// load bitset with flags for function 0x00000007
	if (nIds_ >= 7) {
		f_7_EBX_ = data_[7][1];
		f_7_ECX_ = data_[7][2];
	}
	// Calling __cpuid with 0x80000000 as the function_id argument
	// gets the number of the highest valid extended ID.
	__cpuid(cpui.data(), 0x80000000);
	nExIds_ = cpui[0];
	char brand[0x40];
	memset(brand, 0, sizeof(brand));
	for (int_fast32_t i = 0x80000000; i <= nExIds_; ++i) {
		__cpuidex(cpui.data(), i, 0);
		extdata_.push_back(cpui);
	}
	// load bitset with flags for function 0x80000001
	if (nExIds_ >= 0x80000001) {
		f_81_ECX_ = extdata_[1][2];
		f_81_EDX_ = extdata_[1][3];
	}
	// int_fast32_terpret CPU brand std::string if reported
	if (nExIds_ >= 0x80000004) {
		memcpy(brand, extdata_[2].data(), sizeof(cpui));
		memcpy(brand + 16, extdata_[3].data(), sizeof(cpui));
		memcpy(brand + 32, extdata_[4].data(), sizeof(cpui));
		brand_ = brand;
	}
}

std::string main_CPU() {
	return InstructionSet::Vendor() + " " + InstructionSet::Brand();
}

void additional_CPU_info() {
	auto support_message = [](std::string isa_feature, bool is_supported) {
		std::cout << isa_feature << (is_supported ? " supported" : " not supported") << '\n';
	};
	support_message("3DNOW", InstructionSet::_3DNOW());
	support_message("3DNOWEXT", InstructionSet::_3DNOWEXT());
	support_message("ABM", InstructionSet::ABM());
	support_message("ADX", InstructionSet::ADX());
	support_message("AES", InstructionSet::AES());
	support_message("AVX", InstructionSet::AVX());
	support_message("AVX2", InstructionSet::AVX2());
	support_message("AVX512CD", InstructionSet::AVX512CD());
	support_message("AVX512ER", InstructionSet::AVX512ER());
	support_message("AVX512F", InstructionSet::AVX512F());
	support_message("AVX512PF", InstructionSet::AVX512PF());
	support_message("BMI1", InstructionSet::BMI1());
	support_message("BMI2", InstructionSet::BMI2());
	support_message("CLFSH", InstructionSet::CLFSH());
	support_message("CMPXCHG16B", InstructionSet::CMPXCHG16B());
	support_message("CX8", InstructionSet::CX8());
	support_message("ERMS", InstructionSet::ERMS());
	support_message("F16C", InstructionSet::F16C());
	support_message("FMA", InstructionSet::FMA());
	support_message("FSGSBASE", InstructionSet::FSGSBASE());
	support_message("FXSR", InstructionSet::FXSR());
	support_message("HLE", InstructionSet::HLE());
	support_message("INVPCID", InstructionSet::INVPCID());
	support_message("LAHF", InstructionSet::LAHF());
	support_message("LZCNT", InstructionSet::LZCNT());
	support_message("MMX", InstructionSet::MMX());
	support_message("MMXEXT", InstructionSet::MMXEXT());
	support_message("MONITOR", InstructionSet::MONITOR());
	support_message("MOVBE", InstructionSet::MOVBE());
	support_message("MSR", InstructionSet::MSR());
	support_message("OSXSAVE", InstructionSet::OSXSAVE());
	support_message("PCLMULQDQ", InstructionSet::PCLMULQDQ());
	support_message("POPCNT", InstructionSet::POPCNT());
	support_message("PREFETCHWT1", InstructionSet::PREFETCHWT1());
	support_message("RDRAND", InstructionSet::RDRAND());
	support_message("RDSEED", InstructionSet::RDSEED());
	support_message("RDTSCP", InstructionSet::RDTSCP());
	support_message("RTM", InstructionSet::RTM());
	support_message("SEP", InstructionSet::SEP());
	support_message("SHA", InstructionSet::SHA());
	support_message("SSE", InstructionSet::SSE());
	support_message("SSE2", InstructionSet::SSE2());
	support_message("SSE3", InstructionSet::SSE3());
	support_message("SSE4.1", InstructionSet::SSE41());
	support_message("SSE4.2", InstructionSet::SSE42());
	support_message("SSE4a", InstructionSet::SSE4a());
	support_message("SSSE3", InstructionSet::SSSE3());
	support_message("SYSCALL", InstructionSet::SYSCALL());
	support_message("TBM", InstructionSet::TBM());
	support_message("XOP", InstructionSet::XOP());
	support_message("XSAVE", InstructionSet::XSAVE());
}