#include "av_vel.h"

AvVel::AvVel()
{
    dvx = new double* [parameters::num_thread];

    for (int i = 0; i < parameters::num_thread; i++)
    {
        dvx[i] = new double[parameters::dim];
    }
}

AvVel::~AvVel()
{
    for (int i = 0; i < parameters::num_thread; i++)
    {
        delete[] dvx[i];
    }
    delete[] dvx;
}

double** AvVel::findAvVel(Particle* particles, int n_total,int niac, int* pair_i, int* pair_j, double* w, double** av)
{
    double epsilon = parameters::epsilon;
    //double* dvx = new double[parameters::dim];;

    for (int i = 0; i < n_total; i++)
    {
        for (int d = 0; d < parameters::dim; d++)
        {
            av[d][i] = 0.;
        }
    }

#pragma omp parallel 
    {
#pragma omp for  schedule(static)
        for (int k = 0; k < niac; k++)
        {
            int i = pair_i[k];
            int j = pair_j[k];

            for (int d = 0; d < parameters::dim; d++)
            {
                dvx[omp_get_thread_num()][d] = particles[i].vx[d] - particles[j].vx[d];
                av[d][i] = av[d][i] - 2 * particles[j].mass * dvx[omp_get_thread_num()][d] / (particles[i].rho + particles[j].rho) * w[k];
                av[d][j] = av[d][j] + 2 * particles[i].mass * dvx[omp_get_thread_num()][d] / (particles[i].rho + particles[j].rho) * w[k];
            }
        }
#pragma omp for schedule(static)  
        for (int i  = 0; i < n_total; i++)
        {
            for (int d = 0; d < parameters::dim; d++)
            {
                av[d][i] = epsilon * av[d][i];
            }
        }

    }


   
    return av;
}