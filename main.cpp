#include "sand.h"
#include <stdio.h>
#include <Eigen/Dense>
#include <iostream>

using Eigen::MatrixXd;
typedef Eigen::Matrix<float, 3, 3> matrix;

int
main() {
	float dt = 0.1;
	float h = 0.01;
	Wall zp(0, 0, 1, 0, 1, 0, 0, 1, 0.01);
	Wall xp(1, 0, 0, 0, 1, 0, 0, 1, 0.01);
	Wall yp(0, 1, 0, 0, 1, 0, 0, 1, 0.01);
	Sand sand(10, 10, 10, 10, h, dt);
	sand.random_initial_positions();
	sand.update_weight();
	sand.add_obstacle(zp);
	sand.add_obstacle(xp);
	sand.add_obstacle(yp);
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