#include "sand.h"
#include <stdio.h>
#include <Eigen/Dense>
#include <iostream>

using Eigen::MatrixXd;
typedef Eigen::Matrix<float, 3, 3> matrix;

int
main() {
	Sand sand(100, 100, 100, 100, 100);
	matrix Sigma = matrix::Identity();
	Sigma = Sigma * 2;
	std::cout << Sigma << std::endl;
	float a = 1;
	matrix expH = matrix::Identity();
	float delgam;
	sand.project(Sigma, a, expH, delgam);

	matrix deriv;
	sand.energy_derivative(Sigma, deriv);



	return 0;
}