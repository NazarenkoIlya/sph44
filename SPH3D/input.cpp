#include "input.h"
void ParticleInitialization::LoadShockTube(Particle* particles, int* n_total)
{
    double space_x;
    *n_total = 400;

    space_x = 0.6 / 80;
    for (int i = 0; i < *n_total; i++)
    {
        particles[i].mass = 0.75 / 400;
        particles[i].hsml = 0.015;
        particles[i].type = 1;
        for (int d = 0; d < parameters::dim; d++)
        {
            particles[i].vx[d] = 0.0;
        }
    }
    for (int i = 0; i < 320; i++)
    {
        particles[i].x[0] = -0.6 + space_x / 4 * i;
    }
    for (int i = 320; i < *n_total; i++)
    {
        particles[i].x[0] = space_x * (i - 319);
    }


    for (int i = 0; i < *n_total; i++)
    {
        if (particles[i].x[0] <= 1.0E-8)
        {
            particles[i].u = 2.5;
            particles[i].rho = 1;
            particles[i].p = 1;
        }
        if (particles[i].x[0] > 1.0E-8)
        {
            particles[i].u = 1.795;;
            particles[i].rho = 0.25;
            particles[i].p = 0.1795;
        }
    }
}
void ParticleInitialization::LoadShearCavity(Particle* particles, int* n_total)
{
    double xl, yl, dx, dy;
    int mp = 40, np = 40;
    *n_total = mp * np;
    xl = 1.e-3f;
    yl = 1.e-3;
    dx = xl / mp;
    dy = yl / np;

    for (int i = 0; i < mp; i++)
    {
        for (int j = 0; j < np; j++)
        {
            int k = j + i * np;
            particles[k].x[0] = i * dx + dx / 2;
            particles[k].x[1] = j * dy + dy / 2;
        }
    }

    for (int i = 0; i < mp * np; i++)
    {
        particles[i].vx[0] = 0.;
        particles[i].vx[1] = 0.;
        particles[i].rho = 1000;
        particles[i].mass = dx * dy * particles[i].rho;
        particles[i].p = 0.;
        particles[i].u = 357.1;
        particles[i].type = 2;
        particles[i].hsml = dx;
    }
}
void ParticleInitialization::LoadDambBroke(Particle* particles, int* n_total)
{
    double xl, yl, dx, dy;
    int mp = 100, np = 200;
    *n_total = mp * np;
    xl = 0.5e-3f;
    yl = 1.e-3f;
    dx = xl / mp;
    dy = yl / np;

    parameters::diff = dx;

    for (int i = 0; i < mp; i++)
    {
        for (int j = 0; j < np; j++)
        {
            int k = j + i * np;
            particles[k].x[0] = i * dx + dx / 2.;
            particles[k].x[1] = j * dy + dy / 2.;
        }
    }
    //*n_total = mp * np;
    //for (int i = 0; i < mp; i++)
    //{
    //    for (int j = 0; j < np; j++)
    //    {
    //        int k = j + i * np;
    //        particles[k].x[0] = i * dx + dx / 2.;
    //        particles[k].x[1] = j * dy + dy / 2.;
    //    }
    //}


    for (int i = 0; i < mp * np; i++)
    {
        particles[i].vx[0] = 0.;
        particles[i].vx[1] = 0.;
        particles[i].rho = 1000;
        particles[i].mass = dx * dy * particles[i].rho;
        particles[i].p = 0.;
        particles[i].u = 357.1;
        particles[i].type = 2;
        particles[i].hsml = dx;
    }
}
void ParticleInitialization::spitString(std::string str, std::string strs[])
{
    std::string delimiter = " ";
    size_t pos = 0;
    int n = 0;
    while ((pos = str.find(delimiter)) != std::string::npos)
    {
        strs[n] = str.substr(0, pos);
        str.erase(0, pos + delimiter.length());
        n++;
    }
}
void ParticleInitialization::loadInitialParticle(Particle *particles, int *n_total)
{
    std::string str, strs[100];
    int im;
    if (parameters::confingInput)
    {
        std::ifstream sr1(parameters::_path + "ini_xv.txt");
        std::ifstream sr2(parameters::_path + "ini_state.txt");
        std::ifstream sr3(parameters::_path + "ini_other.txt");
       
        std::getline(sr1, str);
        *n_total = atoi(str.c_str());
        //sr1 >> *n_total;
        std::string delimiter = " ";

        for (int i = 0; i < *n_total; i++)
        {
            //std::string str, strs[100];
            std::getline(sr1, str);
            spitString(str, strs);
            int n = 0;

            im = atoi(strs[n].c_str());
            for (int d = 0; d < parameters::dim; d++)
            {
                particles[i].x[d] = atof(strs[++n].c_str());
            }
            for (int d = 0; d < parameters::dim; d++)
            {
                particles[i].vx[d] = atof(strs[++n].c_str());
            }

            std::getline(sr2, str);
            spitString(str, strs);

            im = atoi(strs[0].c_str());
            particles[i].mass = atof(strs[1].c_str());
            particles[i].rho = atof(strs[2].c_str());
            particles[i].p = atof(strs[3].c_str());
            particles[i].u = atof(strs[4].c_str());

            std::getline(sr3, str);
            spitString(str, strs);

            im = atoi(strs[0].c_str());
            particles[i].type = atoi(strs[1].c_str());
            particles[i].hsml = atof(strs[2].c_str());
        }
        sr1.close();
        sr2.close();
        sr3.close();
    }
    else
    {
        std::ofstream sw1(parameters::_path + "ini_xv.txt");
        std::ofstream sw2(parameters::_path + "ini_state.txt");
        std::ofstream sw3(parameters::_path + "ini_other.txt");

        if (parameters::shockTube) LoadShockTube(particles, n_total);
        if (parameters::shearCavity) LoadShearCavity(particles, n_total);
        if (parameters::dambBrake) LoadDambBroke(particles, n_total);

        sw1 << *n_total << std::endl;
        for (int i = 0; i < *n_total; i++)
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
}
ParticleInitialization::ParticleInitialization(){}
ParticleInitialization::~ParticleInitialization(){}
