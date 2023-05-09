#include "density.h"
#include <iostream>
#include <fstream>
Density::Density()
{
    hv = new double* [parameters::num_thread];

    for (int i = 0; i < parameters::num_thread; i++)
    {
        hv[i] = new double[parameters::dim];
    }
    wi = new double[parameters::maxn];
    dvx = new double[parameters::dim];
}

Density::~Density()
{
    for (int i = 0; i < parameters::num_thread; i++)
    {
        delete[] hv[i];
    }
    delete[] hv;
    //delete[] hv;
    delete[] wi;
    delete[] dvx;
}


void Density::findSumDensity(int n_total, int niac, Particle* particles, int* pair_i, int* pair_j, double* w)
{
    double selfdens = 0.f;
    double r = 0.f;
    Kernel kernel;

#pragma omp parallel
    {
       
#pragma omp  for private(selfdens,kernel) schedule(static)
        for (int i = 0; i < n_total; i++)
        {
            kernel.findKernel(r, hv[omp_get_thread_num()], particles[i].hsml, &selfdens, hv[omp_get_thread_num()]);
            wi[i] = selfdens * particles[i].mass / particles[i].rho;
        }

#pragma omp  for  schedule(static)  
        for (int k = 0; k < niac; k++)
        {
            int i = pair_i[k];
            int j = pair_j[k];

            wi[i] = wi[i] + particles[j].mass / particles[j].rho * w[k];
            wi[j] = wi[j] + particles[i].mass / particles[i].rho * w[k];
        }

#pragma omp  for private(selfdens,kernel) schedule(static)
        for (int i = 0; i < n_total; i++)
        {
            kernel.findKernel(r, hv[omp_get_thread_num()], particles[i].hsml, &selfdens, hv[omp_get_thread_num()]);
            particles[i].rho = selfdens * particles[i].mass;
        }
     
#pragma omp for schedule(static)  
        for (int k = 0; k < niac; k++)
        {
            int i = pair_i[k];
            int j = pair_j[k];

            particles[i].rho += particles[j].mass * w[k];
            particles[j].rho += particles[i].mass * w[k];
        }
        if (parameters::norDensity)
        {
#pragma  omp for schedule(static)  
            for (int i = 0; i < n_total; i++)
            {
                particles[i].rho /= wi[i];
            }
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