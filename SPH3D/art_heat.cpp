#include"art_heat.h"
ArtHeat::ArtHeat()
{
    dvx = new double[parameters::dim];
    vcc = new double[parameters::maxn];
    dedt = new double[parameters::maxn];

    g1 = 0.1;
    g2 = 1.0;
 
}

ArtHeat::~ArtHeat()
{
    delete[] vcc;
    delete[] dvx;
    delete[] dedt;
}

double* ArtHeat::find(int n_total,int niac, int* pair_i, int* pair_j, double* w, double** dwdx, Particle* particles)
{
    double dx, rr, h, mrho, mhsml, hvcc, mui, muj, muij, rdwdx;
#pragma omp parallel for  
    for (int i = 0; i < n_total; i++)
    {
        vcc[i] = 0.f;
        dedt[i] = 0.f;
    }
#pragma omp parallel for private (hvcc)
    for (int k = 0; k < niac; k++)
    {
        int i = pair_i[k];
        int j = pair_j[k];

        for (int d = 0; d < parameters::dim; d++)
        {
            dvx[d] = particles[j].vx[d] - particles[i].vx[d];
        }
        hvcc = dvx[0] * dwdx[0][k];
        for (int d = 1; d < parameters::dim; d++)
        {
            hvcc += dvx[d] * dwdx[d][k];
        }

        vcc[i] += particles[j].mass * hvcc / particles[j].rho;
        vcc[j] += particles[i].mass * hvcc / particles[i].rho;
    }

#pragma omp parallel for private(hvcc,mhsml,rr,rdwdx,dx,mui,muj,muij,h)
    for (int k = 0; k < niac; k++)
    {
        int i = pair_i[k];
        int j = pair_j[k];

        mhsml = 0.5 * (particles[i].hsml + particles[j].hsml);
        mrho = 0.5f * (particles[i].rho + particles[j].rho);

        rr = 0.0f;
        rdwdx = 0.0f;

        for (int d = 0; d < parameters::dim; d++)
        {
            dx = particles[i].x[d] - particles[j].x[d];
            rr = rr + dx * dx;
            rdwdx = rdwdx + dx * dwdx[d][k];
        }

        mui = g1 * particles[i].hsml * particles[i].c + g2 * particles[i].hsml * particles[i].hsml * (fabs(vcc[i]) - vcc[i]);
        muj = g1 * particles[j].hsml * particles[j].c + g2 * particles[j].hsml * particles[j].hsml * (fabs(vcc[j]) - vcc[j]);

        muij = 0.5 * (mui + muj);

        h = muij / (mrho * (rr + 0.01 * mhsml * mhsml)) * rdwdx;
        dedt[i] += particles[j].mass * h * (particles[i].u - particles[j].u);
        dedt[j] += particles[i].mass * h * (particles[j].u - particles[i].u);

    }
#pragma omp parallel for 
    for (int i = 0; i < n_total; i++)
    {
        dedt[i] = 2.f * dedt[i];
    }

    return dedt;
}