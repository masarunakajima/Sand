#include "sand.h"
#include <stdio.h>
#include <stdlib.h>


typedef Eigen::Matrix<float, 3, 3> matrix;
typedef Eigen::Matrix<float, 3, 1> vector;

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

// MN:in progress
Sand::Sand(int xRes, int yRes, int zRes, int nparticle, int ngrid) {
	xres = xRes;
	yres = yRes;
	zres = zRes;
	n_particle = nparticle;
	n_grid = ngrid;
	initialize();
}


// MN: inprogress
void
Sand::initialize() {

	//float* mp, * Vp0, * ap, apn, * qp, * qpn;
	mp = (float*)malloc(n_particle * sizeof(float));
	Vp0 = (float*)malloc(n_particle * sizeof(float));
	ap = (float*)malloc(n_particle * sizeof(float));
	apn = (float*)malloc(n_particle * sizeof(float));
	qp = (float*)malloc(n_particle * sizeof(float));
	qpn = (float*)malloc(n_particle * sizeof(float));
	//Eigen::Matrix<float, 3, 3>* Bp, * Bpn, * Fp, * Fpn, * FEp, * FEpn, * FPp, * FPpn,
	//	* FEpt, * FPpt, * Cp, * Dp, * Fph, * FEph, * FPph, * gradv;
	Bp = (matrix*)malloc(n_particle * sizeof(matrix));
	Bpn = (matrix*)malloc(n_particle * sizeof(matrix));
	Fp = (matrix*)malloc(n_particle * sizeof(matrix));
	Fpn = (matrix*)malloc(n_particle * sizeof(matrix));
	FEp = (matrix*)malloc(n_particle * sizeof(matrix));
	FEpn = (matrix*)malloc(n_particle * sizeof(matrix));
	FPp = (matrix*)malloc(n_particle * sizeof(matrix));
	FPp = (matrix*)malloc(n_particle * sizeof(matrix));
	FPpn = (matrix*)malloc(n_particle * sizeof(matrix));
	FEpt = (matrix*)malloc(n_particle * sizeof(matrix));
	FPpt = (matrix*)malloc(n_particle * sizeof(matrix));
	Cp = (matrix*)malloc(n_particle * sizeof(matrix));
	Dp = (matrix*)malloc(n_particle * sizeof(matrix));
	Fph = (matrix*)malloc(n_particle * sizeof(matrix));
	FEph = (matrix*)malloc(n_particle * sizeof(matrix));
	FPph = (matrix*)malloc(n_particle * sizeof(matrix));
	gradv = (matrix*)malloc(n_particle * sizeof(matrix));
	Zp = (matrix*)malloc(n_particle * sizeof(matrix));

	//Eigen::Matrix<float, 3, 1>* vp, * vpn, *xp, *xpn, *vpb;
	vp = (vector*)malloc(n_particle * sizeof(vector));
	vpn = (vector*)malloc(n_particle * sizeof(vector));
	xp = (vector*)malloc(n_particle * sizeof(vector));
	xpn = (vector*)malloc(n_particle * sizeof(vector));
	vpb = (vector*)malloc(n_particle * sizeof(vector));
	
	w = (float*)malloc(n_particle * n_grid * sizeof(float));
	gradw = (vector*)malloc(n_particle * n_grid*sizeof(vector));

	xg = (vector*)malloc(n_grid * sizeof(vector));
	vgb = (vector*)malloc(n_grid * sizeof(vector));
	for (int p = 0; p < n_particle; p++) {

		mp[p] = 0.0;
		Vp0[p] = 0.0;
		ap[p] = 0.0;
		apn[p] = 0.0;
		qp[p] = 0.0;
		qpn[p] = 0.0;


		Bp[p] = matrix().Zero();
		Bpn[p] = matrix().Zero();
		Fp[p] = matrix().Zero();
		Fpn[p] = matrix().Zero();
		FEp[p] = matrix().Zero();
		FEpn[p] = matrix().Zero();
		FPp[p] = matrix().Zero();
		FPp[p] = matrix().Zero();
		FPpn[p] = matrix().Zero();
		FEpt[p] = matrix().Zero();
		FPpt[p] = matrix().Zero();
		Cp[p] = matrix().Zero();
		Dp[p] = matrix().Zero();
		Fph[p] = matrix().Zero();
		FEph[p] = matrix().Zero();
		FPph[p] = matrix().Zero();
		gradv[p] = matrix().Zero();
		Zp[p] = matrix().Zero();

		vp[p] = vector().Zero();
		vpn[p] = vector().Zero();
		xp[p] = vector().Zero();
		xpn[p] = vector().Zero();
		vpb[p] = vector().Zero();

	}
	for (int p = 0; p < n_particle; p++) {
		for (int i = 0; i < n_grid; i++) {
			gradw[p * n_particle + i] = vector().Zero();
			w[p * n_particle + i] = 0.0;
		}
	}

	for (int i = 0; i < n_grid; i++) {
		xg[i] = vector().Zero();
		vgb[i] = vector().Zero();
	}

}


//MN: inprogress
int Sand::transer_to_grid() {
	return EXIT_SUCCESS;
}

//MN: complete
int
Sand::update_particle_state() {
	for (int p = 0; p < n_particle; p++) {
		xpn[p] = vector().Zero();
		matrix T = matrix().Zero();
		for (int i = 0; i < n_grid; i++) {
			xpn[p] += w[p*n_particle + i] * 
				(xg[i] + dt * w[p * n_particle + i]*vgb[i]);
			T += vgb[i] * gradw[p * n_particle + i].transpose();
		}
		FEph[p] = (matrix().Identity() + dt * T) * FEp[p];
		FPph[p] = FPp[p];
	}
	return EXIT_SUCCESS;
}


//MN; complete
int
Sand::plasticity_hardening() {
	for (int p = 0; p < n_particle; p++) {
		Eigen::JacobiSVD<matrix> svd = FEpn[p].jacobiSvd();
		matrix U = svd.matrixU();
		matrix V = svd.matrixV();
		matrix Sigma = matrix::Zero();
		for (int i = 0; i < 3; i++)
			Sigma(i, i) = svd.singularValues()[i];
		matrix T;
		float delgam;
		project(Sigma, ap[p], T, delgam);
		FEpn[p] = U * T * V.transpose();
		FPpn[p] = V * T.inverse() * Sigma * V.transpose() * FPph[p];
		qpn[p] = qp[p] + delgam;
		float phi = h0 + (h1 * qpn[p] - h3) * exp(-h2*qpn[p]);
		apn[p] = pow(2 / 3, 0.5) * 2 * sin(phi) / (3 - sin(phi));

	}
	return EXIT_SUCCESS;
}

// MN: inprogress
int
Sand::friction(Eigen::Matrix<float, 3, 1>* v, Eigen::Matrix<float, 3, 1>* gradv) {
	
	return EXIT_SUCCESS;
}

//:MN complete
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




//MN: complete
int
Sand::energy_derivative(const matrix& Sigma, 
	matrix& deriv) {
	matrix inv = Sigma.inverse();
	deriv = inv * Sigma.log() * 2 * mu
		+ lambda * Sigma.log().trace() * inv;
	return EXIT_SUCCESS;
}


// MN:complete
int Sand::force_increment(matrix* F, float b, vector* f) {
	matrix* Ap = (matrix*)malloc(n_particle * sizeof(matrix));
	for (int p = 0; p < n_particle; p++) {
		Eigen::JacobiSVD<matrix> svd = F[p].jacobiSvd();
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
		f[i] = vector::Zero();
		for (int p = 0; p < n_particle; p++) {
			f[i] -= dt / mg[i] * (Ap[p]*gradw[p*n_particle+i]);
		}
		force_increment(Ap, 1, f);
	}
	
	return EXIT_SUCCESS;
}


//MN: complete
int Sand::force_increment_vel(vector* v, vector* f) {
	matrix* Ap = (matrix*)malloc(n_particle * sizeof(matrix));
	for (int p = 0; p < n_particle; p++) {
		matrix T = matrix().Zero();
		for (int i = 0; i < n_grid; i++) {
			T += v[i] * gradw[p * n_particle + i].transpose();
		}
		Ap[p] = (matrix().Identity() + dt * T) * FEp[p];
	}
	return EXIT_SUCCESS;
}
