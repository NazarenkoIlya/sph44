#pragma once
#include <string>
#include "particles.h"
#include "param.h"
#include <fstream>
class ParticleInitialization
{
public:
	void loadInitialParticle(Particle* particles, int* n_total);
	ParticleInitialization();
	~ParticleInitialization();

private:
	void LoadDambBroke(Particle* particles, int* n_total);
	void LoadShockTube(Particle* particles, int* n_total);
	void LoadShearCavity(Particle* particles, int* n_total);
	void spitString(std::string str, std::string strs[]);
};



