#pragma once
#include "param.h"
#include "particles.h"

class ExtForce
{
public:
    ExtForce();
    ~ExtForce();
    double** findExtForce(int nTotal, Particle *particles, int niac, int* pair_i, int* pair_j, double** dvxdt);

private:
   // double* dx;
    double rr0, dd, p1, p2;
    double** dx;
};