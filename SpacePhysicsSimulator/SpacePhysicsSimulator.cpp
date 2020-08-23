#include "C_main.h"
#include <exception>

int_fast32_t main() {
	//try {
		do {} while (C_Main().run());
	//}
	//catch (std::exception t) {
	//	std::cout << t.what() << '\n';
	//	system("pause");
	//	throw;
	//}
	return 0;
}