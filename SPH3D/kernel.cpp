#include "kernel.h"

void Kernel::findKernel(double r, double* dx, double hsml, double* w, double* dwdx)
{
    q = r / hsml;
    *w = 0.f;

    if (parameters::dim == 1)
        factor = 1.f / hsml;
    else if (parameters::dim == 2)
        factor = 15.f / (7.f * M_PI * hsml * hsml);
    else if (parameters::dim == 3)
        factor = 3.f / (2.f * M_PI * hsml * hsml * hsml);

    for (int d = 0; d < parameters::dim; d++)  dwdx[d] = 0.f;
    if (parameters::skf == 1) findCubicSplineKenel(r, dx, hsml, w, dwdx);
    else if (parameters::skf == 2) findGaussKernel(r, dx, hsml, w, dwdx);
    else if (parameters::skf == 3) findQuinticKernel(r, dx, hsml, w, dwdx);
}

Kernel::Kernel()
{
    q = 0., factor = 0.;
}

Kernel::~Kernel(){}

void Kernel::findCubicSplineKenel(double r, double* dx, double hsml, double* w, double* dwdx)
{
    q = r / hsml;
    *w = 0.f;
  
    if (q <= 1.f)
    {
        *w = factor * (2.f / 3.f - (q * q) + (q * q * q) * 0.5f);
        for (int d = 0; d < parameters::dim; d++)
        {
            dwdx[d] = dx[d] * (factor * (-2.f + 3.f * 0.5f * q) / (hsml * hsml));
        }
    }
    else if ( q <= 2)
    {
        *w = factor * (1.f / 6.f * ((2.f - q) * (2.f - q) * (2.f - q)));
        for (int d = 0; d < parameters::dim; d++)
        {
           // dwdx = -diff / dist * (factor * sqr(2.f - q) * 0.5f / hsml);
            dwdx[d] = -(dx[d] / r) *(factor * ((2. - q) * (2. - q)) * 0.5f / hsml);
        }
    }
    else
    {
        *w = 0.f;
        for (int d = 0; d < parameters::dim; d++)
        {
            dwdx[d] = 0.f;
        }
    }
}

inline void Kernel::findGaussKernel(double r, double* dx, double hsml, double* w, double* dwdx)
{
    factor = 1.f / (pow(hsml, parameters::dim) * pow(M_PI, (parameters::dim / 2.)));
    if (q >= 0.0 && q <= 3.0)
    {
        *w = factor * exp(-q * q);
        for (int d = 0; d < parameters::dim; d++)
        {
            dwdx[d] = *w * (-2. * dx[d] / hsml / hsml);
        }
    }
    else
    {
        *w = 0.;
        for (int d = 0; d < parameters::dim; d++)
        {
            dwdx[d] = 0.;
        }
    }
}

inline void Kernel::findQuinticKernel(double r, double* dx, double hsml, double* w, double* dwdx)
{
    if (parameters::dim == 1) factor = 1.f / (120.f * hsml);
    else if (parameters::dim == 2) factor = 7.f / (478.0 * M_PI * hsml * hsml);
    else if (parameters::dim == 3) factor = 1.f / (120.f * M_PI * hsml * hsml * hsml);

    if (q >= 0 && q <= 1)
    {
        *w = factor * ((pow(3 - q, 5) - 6 * pow(2 - q, 5) + 15 * pow(1 - q, 5)));
        for (int d = 0; d < parameters::dim; d++)
        {
            dwdx[d] = factor * ((-120 + 120 * q - 50 * q * q) / (hsml * hsml) * dx[d]);
        }
    }
    else if (q > 1 && q <= 2)
    {
        *w = factor * ((pow(3 - q, 5) - 6 * pow(2 - q, 5)));
        for (int d = 0; d < parameters::dim; d++)
        {
            dwdx[d] = factor * (-5 * pow(3 - q, 4) + 30 * pow(2 - q, 4)) / hsml * (dx[d] / r);
        }
    }
    else if (q > 2 && q <= 3)
    {
        *w = factor * pow(3 - q, 5);
        for (int d = 0; d < parameters::dim; d++)
        {
            dwdx[d] = factor * (-5 * pow(3 - q, 4)) / hsml * (dx[d] / r);
        }
    }
    else
    {
        *w = 0.;
        for (int d = 0; d < parameters::dim; d++)
        {
            dwdx[d] = 0.;
        }
    }
}