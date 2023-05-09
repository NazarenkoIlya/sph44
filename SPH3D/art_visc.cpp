#include "art_visc.h"
#include <iostream>
ArtVisc::ArtVisc()
{
    dvx = new double*[parameters::num_thread];

    for (int i = 0; i < parameters::num_thread; i++)
    {
        dvx[i] = new double[parameters::dim];
    }
    etq = parameters::etq;
    alpha = parameters::alpha;
    beta = parameters::beta;

}

ArtVisc::~ArtVisc()
{
    for (int i = 0; i < parameters::num_thread; i++)
    {
        delete[] dvx[i];
    }
    delete[] dvx;
}

void ArtVisc::find(int n_total, Particle* particles, int niac, int* pair_i, int* pair_j, double** dwdx, double** dvxdt, double* dedt)
{
    etq = parameters::etq;
    alpha = parameters::alpha;
    beta = parameters::beta;
    double dx = 0.f, piv, muv = 0.f, vr = 0.f, rr = 0.f, h = 0.f, mc = 0.f, mrho = 0.f, mhsml = 0.f;

    for (int i = 0; i < n_total; i++)
    {
        for (int d = 0; d < parameters::dim; d++)
        {
            dvxdt[d][i] = 0.f;
        }
        dedt[i] = 0.f;
    }

#pragma omp parallel
    {
#pragma omp  for private(mhsml,vr,rr,dx,muv,mc,h,mrho,piv) schedule(static)
        for (int k = 0; k < niac; k++)
        {
            int i = pair_i[k];
            int j = pair_j[k];

            mhsml = (particles[i].hsml + particles[j].hsml) / 2.;

            vr = 0.f;
            rr = 0.f;

            for (int d = 0; d < parameters::dim; d++)
            {
                dvx[omp_get_thread_num()][d] = particles[i].vx[d] - particles[j].vx[d];
                dx = particles[i].x[d] - particles[j].x[d];
                vr += dvx[omp_get_thread_num()][d] * dx;
                rr += dx * dx;
            }

            if (vr < 0)
            {
                muv = mhsml * vr / (rr + mhsml * mhsml * etq * etq);
                mc = 0.5f * (particles[i].c + particles[j].c);
                mrho = 0.5f * (particles[i].rho + particles[j].rho);
                piv = (beta * muv - alpha * mc) * muv / mrho;

                for (int d = 0; d < parameters::dim; d++)
                {
                    h = -piv * dwdx[d][k];

                    dvxdt[d][i] = dvxdt[d][i] + particles[j].mass * h;
                    dvxdt[d][j] = dvxdt[d][j] - particles[i].mass * h;

                    dedt[i] = dedt[i] - particles[j].mass * dvx[omp_get_thread_num()][d] * h;
                    dedt[j] = dedt[j] - particles[i].mass * dvx[omp_get_thread_num()][d] * h;
                }
            }
        }

#pragma omp  for schedule(static)  
        for (int i = 0; i < n_total; i++)
        {
            dedt[i] = 0.5f * dedt[i];
        }
    }
}