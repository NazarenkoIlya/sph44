#pragma once
#include "param.h"
#include "particles.h"

class HSML
{
public:
	HSML();
	~HSML();
    void upgradeHsml(double dt, int nTotal,Particle * particles, int niac, int* pair_i, int* pair_j, double** dwdx);

private: 
    double fac;
    double* dvx;
    double* vcc;
    double* dhsml;
    void simpleUpgrade(int n_total, Particle* paticles);
    void BenzUpgrade(double dt, int n_total, Particle* particles, int niac, int* pair_i, int* pair_j, double** dwdx);
};
