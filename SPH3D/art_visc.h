#pragma once
#include"particles.h"
#include "param.h"
class ArtVisc
{
public:
	ArtVisc();
	~ArtVisc();
   void find(int n_total,Particle* particles, int niac, int* pair_i, int* pair_j, double** dwdx, double** dvxdt, double* dedt);

private:
    double etq, alpha, beta,  * dvx;
};

