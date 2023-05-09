#pragma once
#include <fstream>
#include "param.h"
class Pressure
{
public:
	Pressure();
	~Pressure();
    void p_art_water(double rho, double* p, double* c);
    void p_gas(double rho, double u, double* p, double* c);

private:

};

