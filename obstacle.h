#pragma once
#ifndef OBSTACLE_
#define OBSTACLE_
#include <iostream>
#include <unsupported/Eigen/MatrixFunctions>
#include <Eigen>



class Obstacle {
public:
	Obstacle() {};
	virtual float phi(Eigen::Matrix<float, 3, 1> x) { return 0; };
	
};

class Wall : public Obstacle {
public:
	float a, b, c, d;
	float denom;
	Wall(float a_, float b_, float c_, float d_) :
		a(a_), b(b_), c(c_), d(d_) {
		denom = pow(a * a + b * b + c * c, 0.5);
	};

	float phi(Eigen::Matrix<float, 3, 1> x);
};

#endif
