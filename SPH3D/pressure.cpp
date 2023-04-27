#include "pressure.h"
Pressure::Pressure(){}

Pressure::~Pressure(){}
void Pressure::p_art_water(double rho, double* p, double* c)
{
    double gamma = 7.0, rho0 = 1000.e0,c0 = 80.;
   // double  b = 1.013e5f;

    double b = c0 * c0 * rho0 / gamma;
    double hrho = rho / rho0;
    if (rho > rho0)  *p = b * (pow(hrho, gamma) - 1);
    else  *p = 0.;
    *c = c0;
   
    //*c = 1480.;

    //*c = 0.01;
   // *p = *c * *c * rho;
}
void Pressure::p_gas(double rho, double u, double* p, double* c)
{
    double gamma = 1.4;
    *p = (gamma - 1) * rho * u;
    *c = sqrt((gamma - 1) * u);
}