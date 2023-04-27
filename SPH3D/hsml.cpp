#include"hsml.h"

HSML::HSML()
{
    double* dvx = new double[parameters::dim];
    double* vcc = new double[parameters::maxn];
    double* dhsml = new double[parameters::maxn];
}

HSML::~HSML()
{
    delete[] dvx;
    delete[] vcc;
    delete[] dhsml;
}

void HSML::BenzUpgrade(double dt, int n_total, Particle* particles, int niac, int* pair_i, int* pair_j, double** dwdx)
{
    double  hvcc;
#pragma omp parallel for
    for (int i = 0; i < n_total; i++)
    {
        vcc[i] = 0.f;
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
            hvcc = hvcc + dvx[d] * dwdx[d][k];
        }
        vcc[i] = vcc[i] + particles[j].mass * hvcc / particles[j].rho;
        vcc[j] = vcc[j] + particles[i].mass * hvcc / particles[i].rho;
    }
#pragma omp parallel for
    for (int i = 0; i < n_total; i++)
    {
        dhsml[i] = (particles[i].hsml / parameters::dim) * vcc[i];
        particles[i].hsml += dt * dhsml[i];

        if (particles[i].hsml <= 0.)
        {
            particles[i].hsml -= dt * dhsml[i];
        }
    }
}
 void HSML::upgradeHsml(double dt, int n_total, Particle* particles, int niac, int* pair_i, int* pair_j, double** dwdx)
{
    if (parameters::sle == 0) {}
    else if (parameters::sle == 1) simpleUpgrade(n_total, particles);
    else if (parameters::sle == 2)  BenzUpgrade(dt, n_total, particles, niac, pair_i, pair_j, dwdx);

}
 void HSML::simpleUpgrade(int n_total, Particle* paticles)
{
    fac = 2.0;
    for (int i = 0; i < n_total; i++)
    {
        paticles[i].hsml = fac * (pow(paticles[i].mass / paticles[i].rho, 1.0 / parameters::dim));
    }
}