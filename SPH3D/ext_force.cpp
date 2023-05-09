#include "ext_force.h"
#include <iostream>

ExtForce::ExtForce()
{
    dx = new double* [parameters::num_thread];

    for (int i = 0; i < parameters::num_thread; i++)
    {
        dx[i] = new double[parameters::dim];
    }
    //dx = new double[parameters::dim];
    //rr0 = 4.9e-5f;
    //dd = 1.e-5f;
    //rr0 = 1.25e-5
    //rr0 = 1.25e-5f;
    rr0 =parameters::rr0;
    dd = parameters::dd;
    p1 = parameters::p1;
    p2 = parameters::p2;

}

ExtForce::~ExtForce()
{
    for (int i = 0; i < parameters::num_thread; i++)
    {
        delete[] dx[i];
    }
    delete[] dx;
    //delete[] dx;
}

 double** ExtForce::findExtForce(int n_total, Particle* particles, int niac, int* pair_i, int* pair_j, double** dvxdt)
{

     double rr, f;
    for (int i = 0; i < n_total; i++)
    {
        for (int d = 0; d < parameters::dim; d++)
        {
            dvxdt[d][i] = 0.f;
        }
    }
#pragma omp parallel 
    {
        if (parameters::selfGravity)
        {
#pragma omp  for schedule(static)  
            for (int i = 0; i < n_total; i++)
            {
                dvxdt[parameters::dim - 1][i] = -9.8f;
            }
        }

#pragma omp  for private (rr,f) schedule(static) 
        for (int k = 0; k < niac; k++)
        {
            int i = pair_i[k];
            int j = pair_j[k];

            if (particles[i].type > 0 && particles[j].type < 0)
            {
                rr = 0.;
                for (int d = 0; d < parameters::dim; d++)
                {
                    dx[omp_get_thread_num()][d] = particles[i].x[d] - particles[j].x[d];
                    rr += dx[omp_get_thread_num()][d] * dx[omp_get_thread_num()][d];
                }
                rr = sqrt(rr);
                if (rr < rr0)
                {
                    f = (pow(rr0 / rr, p1) - pow(rr0 / rr, p2)) / (rr * rr);
                    for (int d = 0; d < parameters::dim; d++)
                    {
                        dvxdt[d][i] += dx[omp_get_thread_num()][d] * f * dd;
                    }
                }
            }
        }
    }
    

    return dvxdt;
}