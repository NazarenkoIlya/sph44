#include "output.h"
#include <iostream>

SavingParticle::SavingParticle()
{
}

SavingParticle::~SavingParticle()
{
}

void SavingParticle::save(int n_total,Particle * particles,std::string name_file1, std::string name_file2, std::string name_file3, std::string iter)
{
        std::ofstream sw1(name_file1 + iter + ".txt");
        std::ofstream sw2(name_file2 + iter + ".txt");
        std::ofstream sw3(name_file3 + iter + ".txt");
        sw1 << n_total << std::endl;
        for (int i = 0; i < n_total; i++)
        {
            sw1 << i << " ";
            for (int d = 0; d < parameters::dim; d++)
            {
                sw1 << particles[i].x[d] << " ";
            }
            for (int d = 0; d < parameters::dim; d++)
            {
                sw1 << particles[i].vx[d] << " ";
            }
            sw1 << std::endl;

            sw2 << i << " " << particles[i].mass << " " << particles[i].rho << " " << particles[i].p << " " << particles[i].u << std::endl;
            sw3 << i << " " << particles[i].type << " " << particles[i].hsml << std::endl;
        }

        sw1.close();
        sw2.close();
        sw3.close();
    
}