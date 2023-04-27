#include "art_visc.h"
#include <iostream>
ArtVisc::ArtVisc()
{
    dvx = new double[parameters::dim];
    etq = 0.1f;
    alpha = 1.f;
    beta = 1.f;

}

ArtVisc::~ArtVisc()
{
    delete[] dvx;
}

void ArtVisc::find(int n_total, Particle* particles, int niac, int* pair_i, int* pair_j, double** dwdx, double** dvxdt, double* dedt)
{
    etq = 0.1f;
    alpha = 1.f;
    beta = 1.f;
    double dx = 0.f, piv, muv = 0.f, vr = 0.f, rr = 0.f, h = 0.f, mc = 0.f, mrho = 0.f, mhsml = 0.f;
//#pragma omp parallel for
    for (int i = 0; i < n_total; i++)
    {
        for (int d = 0; d < parameters::dim; d++)
        {
            dvxdt[d][i] = 0.f;
           // dvx[d] = 0.f;
        }
        dedt[i] = 0.f;
    }
    for (int d = 0; d < parameters::dim; d++)
    {
        dvx[d] = 0.f;
    }
//#pragma omp parallel for private(mhsml,rr,dx,muv,mc,h,mrho,piv) 
    for (int k = 0; k < niac; k++)
    {
        int i = pair_i[k];
        int j = pair_j[k];

        mhsml = (particles[i].hsml + particles[j].hsml) / 2.;

        vr = 0.f;
        rr = 0.f;

        for (int d = 0; d < parameters::dim; d++)
        {
            dvx[d] = particles[i].vx[d] - particles[j].vx[d];
            dx = particles[i].x[d] - particles[j].x[d];
            vr += dvx[d] * dx;
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

                dedt[i] = dedt[i] - particles[j].mass * dvx[d] * h;
                dedt[j] = dedt[j] - particles[i].mass * dvx[d] * h;
            }
        }

    }
  //  std::cout << std::endl;
//#pragma omp parallel  for
    for (int i = 0; i < n_total; i++)
    {
        dedt[i] = 0.5f * dedt[i];
    }
}