#pragma once
#include <utility>
 
class Particle
{
private:
	
public:
	int type;
	double hsml, c, u, p, mass, rho, * x, * vx;
	Particle(int type, double hsml, double c, double u, double p, double mass, double rho, int dim);
	Particle();
	~Particle();
};