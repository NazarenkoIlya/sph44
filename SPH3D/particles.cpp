#include "particles.h"
#include "param.h"

Particle::Particle(int type, double hsml, double c, double u, double p, double mass, double rho, int dim)
{
	this->type = type;
	this->hsml = hsml;
	this->c = c;
	this->u = u;
	this->p = p;
	this->mass = mass;
	this->rho = rho;
	x = new double[dim];
	vx = new double[dim];
}

Particle::Particle()
{
	this->type = 0;
	this->hsml = 0.E0;
	this->c = 0.E0;
	this->u = 0.E0;
	this->p = 0.E0;
	this->mass = 0.E0;
	this->rho = 0.E0;

	x = new double[parameters::dim];
	vx = new double[parameters::dim];
	for (int d = 0; d < parameters::dim; d++)
	{
		x[d] = 0.;
		vx[d] = 0.;
	}
}

Particle::~Particle()
{
	delete[] x;
	delete[] vx;
}
