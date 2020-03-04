#pragma once
#ifndef _ENABLE_EXTENDED_ALIGNED_STORAGE
#define __ENABLE_EXTENDED_ALIGNED_STORAGE
#endif // !_ENABLE_EXTENDED_ALIGNED_STORAGE
#ifndef C_AVX_Vec
#define C_AVX_Vec
#include <cmath>
#include <immintrin.h>
#include <exception>
#include "C_Vector.h"

inline double avx_hypot_3(__m256d vec) noexcept {
	double temp[4];
	_mm256_storeu_pd(temp, _mm256_mul_pd(vec, vec));
	return sqrt(temp[0] + temp[1] + temp[2]);
}

inline double avx_hypot_3_and_pow_3_reciprocal(__m256d vec) noexcept {
	double temp[4];
	_mm256_storeu_pd(temp, _mm256_mul_pd(vec, vec));
	return pow(temp[0] + temp[1] + temp[2], -1.5);
}

class AVX256pd {
private:
	__m256d
		vec256_pd = _mm256_set1_pd(double());
public:
	AVX256pd();
	explicit AVX256pd(double coord_0, double coord_1, double coord_2 = double(), double coord_3 = double()) noexcept;
	AVX256pd(double num);
	AVX256pd(__m256d&& coords_vec);
	~AVX256pd() = default;
	AVX256pd
		&& operator*(AVX256pd vec) noexcept,
		&& operator/(AVX256pd&& divisor) noexcept,
		&& operator/(AVX256pd& vec) noexcept,
		&& operator-(AVX256pd& vec) noexcept,
		&& operator+(AVX256pd&& vec) noexcept,
		& operator-=(AVX256pd&& vec) noexcept,
		& operator+=(AVX256pd vec) noexcept,
		& operator=(double num) noexcept;
	double
		norm(),
		norm_and_pow();
	double
		operator()(AVX256pd& vec_2);
	bool
		operator==(AVX256pd& vec_2);
	operator C_Vector<double>();
	operator __m256d();
	double& operator[](size_t index);
	friend class C_Object;
	friend class C_Analysis;
	friend class C_Universe;
	friend AVX256pd&& fma(AVX256pd& vec, AVX256pd& vec_mul, AVX256pd& vec_add) noexcept;
};

inline double AVX256pd::norm() {
	return avx_hypot_3(vec256_pd);
}

inline double AVX256pd::norm_and_pow() {
	return avx_hypot_3_and_pow_3_reciprocal(vec256_pd);
}

inline double AVX256pd::operator()(AVX256pd& vec_2) {
	return avx_hypot_3(_mm256_sub_pd(vec_2, vec256_pd));
}

inline bool AVX256pd::operator==(AVX256pd& vec_2) {
	uint_fast64_t res[4];
	_mm256_store_pd((double*)res, _mm256_cmp_pd(vec256_pd, vec_2.vec256_pd, _CMP_EQ_OQ));
	return (res[0] & res[1] & res[2] & res[3]) != 0;
}

inline AVX256pd::operator __m256d() {
	return vec256_pd;
}

inline double& AVX256pd::operator[](size_t index) {
	[[likely]] if (index < 4) {
		return reinterpret_cast<double*>(&vec256_pd)[index];
	}
}

inline AVX256pd::operator C_Vector<double>() {
	double temp[4];
	_mm256_storeu_pd(temp, vec256_pd);
	return C_Vector<double>(temp[0], temp[1], temp[2]);
}

inline AVX256pd::AVX256pd() {
	vec256_pd = _mm256_set1_pd(double());
}

inline AVX256pd::AVX256pd(double coord_0, double coord_1, double coord_2, double coord_3) noexcept :
	vec256_pd{ coord_0, coord_1, coord_2, coord_3 } {}

inline AVX256pd::AVX256pd(double num) {
	vec256_pd = _mm256_set1_pd(num);
}

inline AVX256pd::AVX256pd(__m256d&& coords_vec) {
	vec256_pd = coords_vec;
}

inline AVX256pd&& AVX256pd::operator*(AVX256pd vec) noexcept {
	return _mm256_mul_pd(vec256_pd, vec);
}

inline AVX256pd&& AVX256pd::operator/(AVX256pd& vec) noexcept {
	return _mm256_div_pd(vec256_pd, vec);
}

inline AVX256pd&& fma(AVX256pd& vec, AVX256pd& vec_mul, AVX256pd& vec_add) noexcept {
	return _mm256_fmadd_pd(vec, vec_mul, vec_add);
}

inline AVX256pd&& AVX256pd::operator+(AVX256pd&& vec) noexcept {
	return _mm256_add_pd(vec256_pd, vec);
}

inline AVX256pd&& AVX256pd::operator/(AVX256pd&& vec_divisor) noexcept {
	return _mm256_div_pd(vec256_pd, vec_divisor);
}

inline AVX256pd&& AVX256pd::operator-(AVX256pd& vec) noexcept {
	return _mm256_sub_pd(vec256_pd, vec);
}

inline AVX256pd& AVX256pd::operator-=(AVX256pd&& vec) noexcept {
	vec256_pd = _mm256_sub_pd(vec256_pd, vec);
	return *this;
}

inline AVX256pd& AVX256pd::operator+=(AVX256pd vec) noexcept {
	vec256_pd = _mm256_add_pd(vec256_pd, vec);
	return *this;
}

inline AVX256pd& AVX256pd::operator=(double num) noexcept {
	vec256_pd = _mm256_set1_pd(num);
	return *this;
}

#endif // !C_AVX_Vec