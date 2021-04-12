#pragma once
#ifndef SAND_
#define SAND_

#include "obstacle.h"

#include <iostream>
#include <unsupported/Eigen/MatrixFunctions>
#include <Eigen>

typedef float coord[3];
#define DEFAULT_MASS 1;
#define DEFAULT_VOL 0.1;



class Sand {			



public:
	int xres, yres, zres;
	int xyres;
	float* mg,    * vgs,  * xgb, * fg;


	
	float dt = 0.1;
	int h;

	float* sigma;
	float* rho;
	float* v;
	
	float h0, h1, h2, h3;
	float phi;
	int d = 3;
	float lambda;
	float mu;
	int n_particle;
	int n_grid;


	float* w;
	Eigen::Matrix<float, 3, 1>* xg;
	Eigen::Matrix<float, 3, 1>* gradw;


	// particle related objects
	float* mp, * Vp0, * ap, * apn, * qp, * qpn;
	Eigen::Matrix<float, 3, 3>* Bp, * Bpn, * Fp, * Fpn, * FEp, * FEpn, * FPp, *FPpn,
		*FEpt, *FPpt, *Cp, *Dp, *Fph, *FEph, *FPph, *gradv, *Zp;
	Eigen::Matrix<float, 3, 1>* vp, * vpn, *xp, *xpn, *vpb;

	Eigen::Matrix<float, 3, 1>* vg, * vgb, * vgt;

	Eigen::Matrix<float, 3, 1> g = Eigen::Matrix<float, 3, 1>(0, 0, -9.8);

	Sand(int xRes, int yRes, int zRes, int nparticle, float h_, float dt_);




	int transer_to_grid();
	int explicit_grid_step();
	int implicit_grid_step();
	int transfer_to_particles();
	int update_particle_state();
	int plasticity_hardening();
	int friction(Eigen::Matrix<float, 3, 1>* v, Eigen::Matrix<float, 3, 1>* gradv);
	int grid_collisions(Eigen::Matrix<float, 3, 1>* v, Eigen::Matrix<float, 3, 1>* vh);
	int project(const Eigen::Matrix<float, 3, 3>&Sigma, float a, Eigen::Matrix<float, 3, 3>& expH, float& delgam);
	int energy_derivative(const Eigen::Matrix<float, 3, 3>& Sigma, Eigen::Matrix<float, 3, 3>&  deriv);
	int force_increment(Eigen::Matrix<float, 3, 3> * Fp, float b, Eigen::Matrix<float, 3, 1>* f);
	int force_increment_vel(Eigen::Matrix<float, 3, 1>* v, Eigen::Matrix<float, 3, 1>* f);

	void initialize();

	int x2i(Eigen::Matrix<float, 3, 1> x);
	void i2x(int i, Eigen::Matrix<float, 3, 1>& x);
	void random_initial_positions();

	float Nh(float x);
	float dNh(float x);

	void update_weight();
	void get_B(Eigen::Matrix<float, 3, 3>* B);

	//float* lambda;
	//float* G;
	//Float* dG;

	//unsigned short	xres;
	//unsigned short	yres;
	//GzPixel* pixelbuffer;		/* frame buffer array */
	//char* framebuffer;

	//GzCamera		m_camera;
	//short		    matlevel;	        /* top of stack - current xform */
	//GzMatrix		Ximage[MATLEVELS];	/* stack of xforms (Xsm) */
	//GzMatrix		Xnorm[MATLEVELS];	/* xforms for norms (Xim) */
	//GzMatrix		Xsp;		        /* NDC to screen (pers-to-screen) */
	//GzColor		flatcolor;          /* color state for flat shaded triangles */
	//int			interp_mode;
	//int			numlights;
	//GzLight		lights[MAX_LIGHTS];
	//GzLight		ambientlight;
	//GzColor		Ka, Kd, Ks;
	//float		    spec;		/* specular power */
	//GzTexture		tex_fun;    /* tex_fun(float u, float v, GzColor color) */

	//// Constructors
	//GzRender(int xRes, int yRes);
	//~GzRender();

	//// Function declaration

	//// HW1: Display methods
	//int GzDefault();
	//int GzBeginRender();
	//int GzPut(int i, int j, GzIntensity r, GzIntensity g, GzIntensity b, GzIntensity a, GzDepth z);
	//int GzGet(int i, int j, GzIntensity* r, GzIntensity* g, GzIntensity* b, GzIntensity* a, GzDepth* z);

	//int GzFlushDisplay2File(FILE* outfile);
	//int GzFlushDisplay2FrameBuffer();

	//// HW2: Render methods
	//int GzPutAttribute(int numAttributes, GzToken* nameList, GzPointer* valueList);
	//int GzPutTriangle(int numParts, GzToken* nameList, GzPointer* valueList);

	//// HW3
	//int GzPutCamera(GzCamera camera);
	//int GzPushMatrix(GzMatrix	matrix);
	//int GzPopMatrix();

	//// Extra methods: NOT part of API - just for general assistance */
	//inline int ARRAY(int x, int y) { return (x + y * xres); }	/* simplify fbuf indexing */
	//inline short	ctoi(float color) { return(short)((int)(color * ((1 << 12) - 1))); }		/* convert float color to GzIntensity short */


	//// Object Translation
	//int GzRotXMat(float degree, GzMatrix mat);
	//int GzRotYMat(float degree, GzMatrix mat);
	//int GzRotZMat(float degree, GzMatrix mat);
	//int GzTrxMat(GzCoord translate, GzMatrix mat);
	//int GzScaleMat(GzCoord scale, GzMatrix mat);

	//// hW6
	//float Xoffset;
	//float Yoffset;
	//float weight;

};

#endif
