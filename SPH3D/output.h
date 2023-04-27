#pragma once
#include <string>
#include <fstream>
#include "particles.h"
#include "param.h"
class SavingParticle
{
public:
	SavingParticle();
	~SavingParticle();
	void save(int n_total, Particle* particles, std::string name_file1, std::string name_file2, std::string name_file3, std::string iter);
private:

};

