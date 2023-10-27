#include "Matrix.h"
#include <iostream>

int main() {
	setlocale(LC_ALL, "Russian");

	linalg::Matrix k1(5);
	linalg::Matrix k2 = { {111111,3,1, 1, 33333},{2,1,1,333,333},{1,1,1,333,1}, {1, 22, 1,22,22} };
	linalg::Matrix k4 = { {1, 2 ,3},{2,4,6},{-5,3,-1} };
	linalg::Matrix m = { {1,1,1}, {2,1,2}, {3,1,1} };
	linalg::Matrix result = { {-0.5,0.0,0.5}, {2.0,-1.0,0.0}, {-0.5,1.0,-0.5} };
	k4.gauss_forward();
	std::cout <<invert(m);

};
