#include "obstacle.h"

typedef Eigen::Matrix<float, 3, 3> matrix;
typedef Eigen::Matrix<float, 3, 1> vector;

float Wall::phi(vector x) {

	return (a*x(0) + b*x(1)+ c*x(2) + d)/ denom;
}
