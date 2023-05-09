#pragma once
#include "param.h"
#include "particles.h"
class AvVel
{
public:
	AvVel();
	~AvVel();
	double** findAvVel(Particle* particles, int n_total, int niac, int* pair_i, int* pair_j, double* w, double** av);
private:
	double** dvx;
};
