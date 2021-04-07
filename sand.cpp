#include "sand.h"
#include <stdio.h>
#include <stdlib.h>


typedef Eigen::Matrix<float, 3, 3> matrix;

//void initialize_matrix0(matrix m) {
//	for (int i = 0; i < 9; i++) {
//		m[i] = 0;
//	}
//};
//void initialize_matrix1(matrix m) {
//	for (int i = 0; i < 9; i++) {
//		m[i] = 0;
//	}
//	m[2] = 1;
//	m[5] = 1;
//	m[8] = 1;
//};
//


Sand::Sand(int xRes, int yRes, int zRes, int nparticle, int ngrid) {
	xres = xRes;
	yres = yRes;
	zres = zRes;
	n_particle = nparticle;
	n_grid = ngrid;
}

void
Sand::initialize() {

}

int Sand::transer_to_grid() {
	return EXIT_SUCCESS;
}

int Sand::project(const matrix& Sigma, float a, matrix& expH, float& delgam) {
	Sigma.log();
	matrix e = Sigma;
	e = e.log();
	matrix eh = e - matrix::Identity()*e.trace()/d;
	if ((eh == matrix::Zero()) || (e.trace() > 0)) {
		expH = matrix::Identity();
		delgam = 0;
		return EXIT_SUCCESS;
	}
	float norm = eh.norm();
	delgam = norm+(d*lambda+2*mu)/(2*mu)*e.trace()*a;
	if (delgam <= 0) {
		expH = Sigma;
		delgam = 0;
		return EXIT_SUCCESS;
	}
	expH = e - eh * delgam / norm;
	expH = expH.exp();
	return EXIT_SUCCESS;
}

int
Sand::energy_derivative(const matrix& Sigma, 
	matrix& deriv) {
	matrix inv = Sigma.inverse();
	deriv = inv * Sigma.log() * 2 * mu
		+ lambda * Sigma.log().trace() * inv;
	return EXIT_SUCCESS;
}

int Sand::force_increment(matrix* F, float b, matrix* f) {
	matrix* Ap = (matrix*)malloc(n_particle * sizeof(matrix));
	for (int p = 0; p < n_particle; p++) {
		Eigen::JacobiSVD<matrix> svd = F[0].jacobiSvd();
		float delgam;
		matrix Sigma = matrix::Zero();
		for (int i = 0; i < 3; i++)
			Sigma(i, i) = svd.singularValues()[i];
		if (b != 0) {
			project(Sigma, ap[p], Sigma, delgam);
		}
		matrix T;
		energy_derivative(Sigma, T);
		Ap[p] = Vp0[p] * svd.matrixU() * T *
			svd.matrixV().transpose() * Fp[p].transpose();		
	}
	for (int i = 0; i < n_grid; i++) {
		f[i] = matrix::Zero();
		for (int p = 0; p < n_particle; p++) {
			f[i] -= dt / mg[i] * Ap[p];
		}
	}
	return EXIT_SUCCESS;
}
