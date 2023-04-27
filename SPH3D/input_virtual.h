#pragma once
#include <string>
#include "particles.h"
#include "param.h"
#include <fstream>
class InitVirtParticle
{
public:
	InitVirtParticle();
	~InitVirtParticle();
	void createVirtPart(int itimestep, int* nVirt,int n_total, Particle* particles, std::string fileName1, std::string fileName2, std::string fileName3);
private:
	void createVirtDambBreak(int itimestep,int* nVirt, int n_total,Particle* particles);
	void createVirtShearCavity(int* nVirt, int n_total, Particle* particles);
	void spitString(std::string str, std::string strs[]);
};

