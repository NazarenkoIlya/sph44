#pragma once
#include "param.h"
#include "particles.h"
#include <math.h>
# define M_PI           3.14159265358979323846 
//#define _USE_MATH_DEFINES

class Kernel
{
public:
    void findKernel(double r, double* dx, double hsml, double* w, double* dwdx);
   
    Kernel();
    ~Kernel();

private:
    double q , factor;
    void findCubicSplineKenel(double r, double* dx, double hsml, double* w, double* dwdx);
    void findGaussKernel(double r, double* dx, double hsml, double* w, double* dwdx);
    void findQuinticKernel(double r, double* dx, double hsml, double* w, double* dwdx);
};

