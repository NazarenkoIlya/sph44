#include "density.h"
#include <iostream>
#include <fstream>
Density::Density()
{
    hv = new double[parameters::dim];
    wi = new double[parameters::maxn];
    dvx = new double[parameters::dim];
    
    for (int d = 0; d < parameters::dim; d++)
    {
        hv[d] = 0.f;
    }
}

Density::~Density()
{
    delete[] hv;
    delete[] wi;
    delete[] dvx;
}

void Density::findSumDensity(int n_total, int niac, Particle* particles, int* pair_i, int* pair_j, double* w)
{
    for (int d = 0; d < parameters::dim; d++)
    {
        hv[d] = 0.f;
    }
 /*   for (int  i = 0; i < n_total; i++)
    {
        wi[i] = 0.f;
    }*/
    double selfdens = 0.f;
    double r = 0.f;
     //Kernel kernel;
//#pragma omp parallel for private(r,selfdens)
    for (int i = 0; i < n_total; i++)
    {
        r = 0.;
        Kernel kernel;
        kernel.findKernel(r, hv, particles[i].hsml, &selfdens, hv);
        wi[i] = selfdens * particles[i].mass / particles[i].rho;
        
    }
    //for(int i =0;i < 10; i++)
    //{
    //    std::cout << w[i]<< std::endl;
    //    for (int d = 0; d < parameters::dim; d++)
    //    {
    //        std::cout << hv[d] << std::endl;
    //    }
    //}
   
//#pragma omp parallel  for 
    for (int k = 0; k < niac; k++)
    {
        int i = pair_i[k];
        int j = pair_j[k];

        wi[i] = wi[i] + particles[j].mass / particles[j].rho * w[k];
        wi[j] = wi[j] + particles[i].mass / particles[i].rho * w[k];
    }
//#pragma omp parallel  for private(selfdens)
    for (int i = 0; i < n_total; i++)
    {
        Kernel kernel;
        kernel.findKernel(r, hv, particles[i].hsml, &selfdens, hv);
        particles[i].rho = selfdens * particles[i].mass;
    }
//#pragma omp parallel  for
    for (int k = 0; k < niac; k++)
    {
        int i = pair_i[k];
        int j = pair_j[k];

        particles[i].rho += particles[j].mass * w[k];
        particles[j].rho += particles[i].mass * w[k];
    }

    if (parameters::norDensity)
    {
//#pragma  omp parallel for
        for (int i = 0; i < n_total; i++)
        {
            particles[i].rho /= wi[i];
        }
    }
}
 
double* Density::findConDensity(int nTotal, Particle* particles, int niac, int* pair_i, int* pair_j, double** dwdx, double* drhodt)
{
    double vcc;

    for (int i = 0; i < nTotal; i++)
    {
        drhodt[i] = 0.;
    }
#pragma  omp parallel for private (vcc)
    for (int k = 0; k < niac; k++)
    {
        int i = pair_i[k];
        int j = pair_j[k];

        for (int d = 0; d < parameters::dim; d++)
        {
            dvx[d] = particles[i].vx[d] - particles[j].vx[d];
        }
        vcc = dvx[0] * dwdx[0][k];
        for (int d = 1; d < parameters::dim; d++)
        {
            vcc = vcc + dvx[d] * dwdx[d][k];
        }

        drhodt[i] = drhodt[i] + particles[j].mass * vcc;
        drhodt[j] = drhodt[j] + particles[i].mass * vcc;
    }

    return drhodt;
}