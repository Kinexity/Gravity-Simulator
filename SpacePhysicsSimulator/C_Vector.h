#pragma once
#ifndef C_Vec
#define C_Vec
#include <cmath>

template <typename T>
T hypot__(T d1, T d2, T d3) {
	return sqrt(d1 * d1 + d2 * d2 + d3 * d3);
};

template<typename Var_Type = double>
class C_Vector {
private:
	Var_Type
		x[3] = { Var_Type() };
public:
	C_Vector() = default;
	explicit C_Vector(Var_Type a, Var_Type b, Var_Type c) noexcept;
	explicit C_Vector(Var_Type num);
	~C_Vector() = default;
	C_Vector<Var_Type>
		&&operator*(Var_Type) noexcept,
		&&operator/(Var_Type) noexcept,
		&&operator-(C_Vector<Var_Type>) noexcept,
		&&operator+(C_Vector<Var_Type>) noexcept;
	C_Vector<Var_Type>&
		operator+=(C_Vector<Var_Type> vec) noexcept;
	void
		operator-=(C_Vector<Var_Type> vec) noexcept,
		operator=(Var_Type num) noexcept;
	Var_Type
		norm(),
		norm_and_pow();
	Var_Type&
		operator[](size_t index);
	Var_Type
		operator()(C_Vector<Var_Type>& vec_2);
	friend class C_Object;
	friend class C_Analysis;
	friend class C_Universe;
	friend class AVX256pd;
};

template<typename Var_Type>
inline Var_Type hypot_3_and_pow_3_reciprocal(C_Vector<Var_Type>& vec) {
	return pow(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2], -1.5);
}

template<typename Var_Type>
inline Var_Type C_Vector<Var_Type>::norm() {
	return hypot__(x[0], x[1], x[2]);
}

template<typename Var_Type>
inline Var_Type C_Vector<Var_Type>::norm_and_pow()
{
	return hypot_3_and_pow_3_reciprocal(*this);
}

template<typename Var_Type>
inline Var_Type & C_Vector<Var_Type>::operator[](size_t index) {
	return x[index];
}

template<typename Var_Type>
inline Var_Type C_Vector<Var_Type>::operator()(C_Vector<Var_Type>& vec_2) {
	return sqrt(fma((vec_2.x[0] - x[0]), (vec_2.x[0] - x[0]), fma((vec_2.x[1] - x[1]), (vec_2.x[1] - x[1]), (vec_2.x[2] - x[2]) * (vec_2.x[2] - x[2]))));
}

template<typename Var_Type>
C_Vector<Var_Type>::C_Vector(Var_Type a, Var_Type b, Var_Type c) noexcept :
	x{ a,b,c } {}

template<typename Var_Type>
inline C_Vector<Var_Type>::C_Vector(Var_Type num) {
	x[0] = x[1] = x[2] = num;
}

template<typename Var_Type>
C_Vector<Var_Type>&& C_Vector<Var_Type>::operator*(Var_Type m) noexcept {
	return C_Vector<Var_Type>(x[0] * m, x[1] * m, x[2] * m);
}

template<typename Var_Type>
C_Vector<Var_Type>&& C_Vector<Var_Type>::operator+(C_Vector<Var_Type> b) noexcept {
	return C_Vector<Var_Type>(x[0] + b.x[0], x[1] + b.x[1], x[2] + b.x[2]);
}

template<typename Var_Type>
C_Vector<Var_Type>&& C_Vector<Var_Type>::operator/(Var_Type m) noexcept {
	return C_Vector<Var_Type>((x[0] / m), (x[1] / m), (x[2] / m));
}

template<typename Var_Type>
C_Vector<Var_Type>&& C_Vector<Var_Type>::operator-(C_Vector<Var_Type> b) noexcept {
	return C_Vector<Var_Type>(x[0] - b.x[0], x[1] - b.x[1], x[2] - b.x[2]);
}

template<typename Var_Type>
inline C_Vector<Var_Type>& C_Vector<Var_Type>::operator+=(C_Vector<Var_Type> vec) noexcept {
	x[0] += vec.x[0];
	x[1] += vec.x[1];
	x[2] += vec.x[2];
	return *this;
}

template<typename Var_Type>
inline void C_Vector<Var_Type>::operator-=(C_Vector<Var_Type> vec) noexcept {
	x[0] -= vec.x[0];
	x[1] -= vec.x[1];
	x[2] -= vec.x[2];
}

template<typename Var_Type>
void C_Vector<Var_Type>::operator=(Var_Type num) noexcept {
	x[0] = x[1] = x[2] = num;
}

#endif // !C_Vec