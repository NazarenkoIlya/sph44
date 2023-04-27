#pragma once
#include"param.h"
#include "particles.h"

class ArtHeat
{
public:
    ArtHeat();
    ~ArtHeat();
    double* find(int n_total,int niac, int* pair_i, int* pair_j, double* w, double** dwdx, Particle *paticles);

private:
    double  g1, g2, * dvx, * vcc;
    double* dedt;
};