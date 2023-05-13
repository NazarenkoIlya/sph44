#include "input_virtual.h"

void InitVirtParticle::createVirtDambBreak(int itimestep,int* nVirt, int n_total, Particle* particles)
{
    int im, mp;
    double xl, dx, v_inf;
    *nVirt = 0;
    mp = 500;
    xl = 2.5E-3;
    dx = xl / mp;
    v_inf = 1.E-3;

 /*   for (int i = 0; i < 2 * mp + 1; i++)
    {
        particles[n_total + *nVirt].x[0] = i * dx / 2.0;
        particles[n_total + *nVirt].x[1] = xl/2.;
        particles[n_total + *nVirt].vx[0] = v_inf;
        particles[n_total + *nVirt].vx[1] = 0.0;
        *nVirt = *nVirt + 1;
    }*/
    for (int i = 0; i < 2 * mp + 1; i++)
    {
        particles[n_total + *nVirt].x[0] = i * dx / 2.;
        particles[n_total + *nVirt].x[1] = xl/2;
        particles[n_total + *nVirt].vx[0] = 0.;
        particles[n_total + *nVirt].vx[1] = 0.0;
        *nVirt = *nVirt + 1;
    }
    for (int i = 0; i < 2 * mp + 1; i++)
    {
        particles[n_total + *nVirt].x[0] = i * dx / 2;
        particles[n_total + *nVirt].x[1] = 0.0;
        particles[n_total + *nVirt].vx[0] = 0.;
        particles[n_total + *nVirt].vx[1] = 0.0;
        *nVirt = *nVirt + 1;
    }
    for (int i = 0; i < 2 * mp/2 - 1; i++)
    {
        particles[n_total + *nVirt].x[0] = 0.0;
        particles[n_total + *nVirt].x[1] = (i + 1) * dx / 2;
        particles[n_total + *nVirt].vx[0] = 0.0;
        particles[n_total + *nVirt].vx[1] = 0.0;
        *nVirt = *nVirt + 1;
    }
    //if (itimestep  == 0)
    //{
    //    for (int i = 0; i < 2 * mp - 1; i++)
    //    {
    //        particles[n_total + *nVirt].x[0] = xl/2;
    //        particles[n_total + *nVirt].x[1] = (i + 1) * dx / 2;
    //        particles[n_total + *nVirt].vx[0] = 0.0;
    //        particles[n_total + *nVirt].vx[1] = 0.0;
    //        *nVirt = *nVirt + 1;
    //    }
    //}
    for (int i = mp/10; i < 2 * mp / 2 - 1; i++)
    {
        particles[n_total + *nVirt].x[0] = 0.51e-3;
        particles[n_total + *nVirt].x[1] = (i + 1) * dx / 2;
        particles[n_total + *nVirt].vx[0] = 0.0;
        particles[n_total + *nVirt].vx[1] = 0.0;
        *nVirt = *nVirt + 1;
    }
    for (int i = 0; i < 2 * mp/2 - 1; i++)
    {
        particles[n_total + *nVirt].x[0] = xl;
        particles[n_total + *nVirt].x[1] = (i + 1) * dx / 2;
        particles[n_total + *nVirt].vx[0] = 0.0;
        particles[n_total + *nVirt].vx[1] = 0.0;
        *nVirt = *nVirt + 1;
    }
    for (int i = 0; i < *nVirt; i++)
    {
        particles[n_total + i].rho = 1000;
        particles[n_total + i].mass = particles[n_total + i].rho * dx * dx;
        particles[n_total + i].p = 0.0;
        particles[n_total + i].u = 357.1;
        particles[n_total + i].type = -2;
        particles[n_total + i].hsml = dx;
    }
}
void InitVirtParticle::createVirtShearCavity(int* nVirt, int n_total, Particle* particles)
{
    int im, mp;
    double xl, dx, v_inf;
    *nVirt = 0;
    mp = 40;
    xl = 1E-3;
    dx = xl / mp;
    v_inf = 1E-3;

    for (int i = 0; i < 2 * mp + 1; i++)
    {
        particles[n_total + *nVirt].x[0] = i * dx / 2;
        particles[n_total + *nVirt].x[1] = xl;
        particles[n_total + *nVirt].vx[0] = v_inf;
        particles[n_total + *nVirt].vx[1] = 0.0;
        *nVirt = *nVirt + 1;
    }
    for (int i = 0; i < 2 * mp + 1; i++)
    {
        particles[n_total + *nVirt].x[0] = i * dx / 2;
        particles[n_total + *nVirt].x[1] = 0.0;
        particles[n_total + *nVirt].vx[0] = 0.0;
        particles[n_total + *nVirt].vx[1] = 0.0;
        *nVirt = *nVirt + 1;
    }
    for (int i = 0; i < 2 * mp - 1; i++)
    {
        particles[n_total + *nVirt].x[0] = 0.0;
        particles[n_total + *nVirt].x[1] = (i + 1) * dx / 2;
        particles[n_total + *nVirt].vx[0] = 0.0;
        particles[n_total + *nVirt].vx[1] = 0.0;
        *nVirt = *nVirt + 1;
    }
    for (int i = 0; i < 2 * mp - 1; i++)
    {
        particles[n_total + *nVirt].x[0] = xl;
        particles[n_total + *nVirt].x[1] = (i + 1) * dx / 2;
        particles[n_total + *nVirt].vx[0] = 0.0;
        particles[n_total + *nVirt].vx[1] = 0.0;
        *nVirt = *nVirt + 1;
    }
    for (int i = 0; i < *nVirt; i++)
    {
        particles[n_total + i].rho = 1000;
        particles[n_total + i].mass = particles[n_total + i].rho * dx * dx;
        particles[n_total + i].p = 0.0;
        particles[n_total + i].u = 357.1;
        particles[n_total + i].type = -2;
        particles[n_total + i].hsml = dx;
    }
}
void InitVirtParticle::spitString(std::string str, std::string strs[])
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
void InitVirtParticle::createVirtPart(int itimestep, int* nVirt, int n_total, Particle* particles, std::string fileName1, std::string fileName2, std::string fileName3)
{/*std::ifstream sr1("xv_vp.txt");
std::ifstream sr2("state_vp.txt");
std::ifstream sr3("other_vp.txt");*/
    int im;
    if (parameters::vpInput)
    {
        std::string str, strs[100];
        std::ifstream sr1(fileName1);
        std::ifstream sr2(fileName2);
        std::ifstream sr3(fileName3);
       
        std::getline(sr1, str);
        *nVirt = atoi(str.c_str());
        //sr1 >> *nVirt;

        for (int j = 0; j < *nVirt; j++)
        {
            int i = n_total + j;
           // std::string str, strs[100];
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
        if (parameters::dambBrake) createVirtDambBreak(itimestep,nVirt,n_total,particles);
        if (parameters::shearCavity) createVirtShearCavity(nVirt, n_total, particles);
    }
    if ((itimestep % parameters::saveStep) == 0)
    {
        std::ofstream sw1(fileName1);
        std::ofstream sw2(fileName2);
        std::ofstream sw3(fileName3);
        sw1 << *nVirt << std::endl;
        for (int i = n_total; i < n_total + *nVirt; i++)
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

            sw2 << i << " " << particles[i].mass << " " << particles[i].rho << " " << particles[i].p << " " << particles[i].u << " " << std::endl;
            sw3 << i << " " << particles[i].type << " " << particles[i].hsml << " " << std::endl;
        }

        sw1.close();
        sw2.close();
        sw3.close();
    }
}
InitVirtParticle::InitVirtParticle() {}
InitVirtParticle::~InitVirtParticle() {}