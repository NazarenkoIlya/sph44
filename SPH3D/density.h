#pragma once
#include "param.h"
#include "particles.h"
#include "kernel.h"
class Density
{
public:
	Density();
	~Density();
    void findSumDensity(int n_total, int niac,Particle *particles, int* pair_i, int* pair_j, double* w);
    double* findConDensity(int nTotal, Particle* particles, int niac, int* pair_i, int* pair_j, double** dwdx, double* drhodt);


private:
    
    double  * hv, * wi;
    double* dvx;
   
};

