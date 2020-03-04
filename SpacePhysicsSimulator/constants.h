#pragma once
#ifndef consts
#define consts
#include <iostream>
#include <string>
#include <memory>
#include <utility>
#include <functional>
#include <limits>
#include "C_Vector.h"
#include "PureCPPLib/io.h"
#include "PureCPPLib/C_Event_Log.h"
#include "C_Universe.h"
#include "C_Sim_Basic_Data.h"
using namespace std;
#undef max
/*
void check_for_cuda_error(cudaError_t state, PCL::C_Event_Log& ev_log_obj, int_fast32_t line_number, string file);

#define CHECK_FOR_CUDA_ERROR(fn) check_for_cuda_error(fn, event_log_obj, __LINE__, __FILE__)

const dim3
block_dim = { 16, 1, 1 };
const dim3
warp_dim = { 16, 1, 1 };
const uint_fast64_t
num_of_blocks = static_cast<uint_fast64_t>(warp_dim.x) * static_cast<uint_fast64_t>(warp_dim.y) * static_cast<uint_fast64_t>(warp_dim.z);
const uint_fast64_t
size_of_block = static_cast<uint_fast64_t>(block_dim.x) * static_cast<uint_fast64_t>(block_dim.y) * static_cast<uint_fast64_t>(block_dim.z);
constexpr uint_fast64_t
num_of_days_to_analyse = 1;
const uint_fast64_t
grid_size = size_of_block * num_of_blocks;*/
constexpr uint_fast64_t
num_of_analysis_types = 2;
#define nm(x) #x << " = " << x
#endif // !consts