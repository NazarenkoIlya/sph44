#pragma once
#include"param.h"
#include "particles.h"
#include "pressure.h"

class IntForce
{
public:
	IntForce();
	~IntForce();
    void find(int nTotal, Particle * particles, int niac, int* pair_i, int* pair_j, double** dwdx, double** dvxdt, double* tdsdt, double* dedt, double* eta);

private:
    double** dvx;
   // double* dvx;
    double* txx;
    double* tyy;
    double* txy;
    void init(int n_total, double* dedt, double** dvxdt, double* tdsdt);
};
