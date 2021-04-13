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
	virtual void gradphi(Eigen::Matrix<float, 3, 1> x, Eigen::Matrix<float, 3, 1> gp) 
	{ };

	Obstacle(bool sticky_, bool separating_, bool slipping_, float mu_) :
		sticky(sticky_), separating(separating_), slipping(slipping_),
	mu(mu_){};

	bool is_sticky() { return sticky; }
	bool is_separating() { return separating; }
	bool is_slipping() { return slipping; }
	float get_mu() { return mu; }
private:
	bool sticky, separating, slipping;
	float mu;

};

class Wall : public Obstacle {
public:

	
	Wall(float a_, float b_, float c_, float d_, float coef_,
		bool sticky_, bool separating_, bool slipping_, float mu_) :
		Obstacle(sticky_, separating_, slipping_, mu_),
		a(a_), b(b_), c(c_), d(d_), coef(coef_)
	{
		
		denom = pow(a * a + b * b + c * c, 0.5);
		gp(0) = coef * a / denom;
		gp(1) = coef * b / denom;
		gp(2) = coef * c / denom;
	};
	
	float phi(Eigen::Matrix<float, 3, 1> x);
	void gradphi(Eigen::Matrix<float, 3, 1> x, Eigen::Matrix<float, 3, 1>& gp);

private:
	float coef;
	float a, b, c, d;
	float denom;
	Eigen::Matrix<float, 3, 1> gp;
};

#endif
