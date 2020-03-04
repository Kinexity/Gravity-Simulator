#pragma once
#ifndef C_Bs
#define C_Bs

class C_Base {
public:
	C_Base() = default;
	virtual ~C_Base() = default;
	virtual void run() = 0;
};

#endif // !C_Bs