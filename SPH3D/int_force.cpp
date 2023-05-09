#include"int_force.h"

IntForce::IntForce()
{
    dvx = new double* [parameters::num_thread];

    for (int i = 0; i < parameters::num_thread; i++)
    {
        dvx[i] = new double[parameters::dim];
    }
   // dvx = new double[parameters::dim];
    txx = new double[parameters::maxn];
    tyy = new double[parameters::maxn];
    txy = new double[parameters::maxn];

}

IntForce::~IntForce()
{
    for (int i = 0; i < parameters::num_thread; i++)
    {
        delete[] dvx[i];
    }
    delete[] dvx;
  //  delete[] dvx;
    delete[] txx;
    delete[] tyy;
    delete[] txy;
}
void IntForce::init(int n_total, double* dedt, double** dvxdt,double *tdsdt)
{
    for (int i = 0; i < n_total; i++)
    {
        txx[i] = 0.f;
        tyy[i] = 0.f;
        txy[i] = 0.f;
        tdsdt[i] = 0.f;
        //vcc[i] = 0.f;
        dedt[i] = 0.f;
        for (int d = 0; d < parameters::dim; d++)
        {
            dvxdt[d][i] = 0.f;
            //std::cout << dvxdt[d][i] << std::endl;
        }
    }
}


void IntForce::find(int n_total, Particle* particles, int niac, int* pair_i, int* pair_j, double** dwdx, double** dvxdt, double* tdsdt, double* dedt, double* eta)
{
    init(n_total, dedt, dvxdt,tdsdt);
    double hxx, hyy, hxy, h, he, rhoij;
    Pressure pressure;
#pragma omp parallel 
    {
        if (parameters::visc)
        {
#pragma omp for private(hxx, hyy, hxy) schedule(static)  
            for (int k = 0; k < niac; k++)
            {
                int i = pair_i[k];
                int j = pair_j[k];

                for (int d = 0; d < parameters::dim; d++)
                {
                    dvx[omp_get_thread_num()][d] = particles[j].vx[d] - particles[i].vx[d];
                }
                if (parameters::dim == 1)
                {
                    hxx = 2.f * dvx[omp_get_thread_num()][0] * dwdx[0][k];
                }
                else if (parameters::dim == 2)
                {
                    hxx = 2.f * dvx[omp_get_thread_num()][0] * dwdx[0][k] - dvx[omp_get_thread_num()][1] * dwdx[1][k];
                    hxy = dvx[omp_get_thread_num()][0] * dwdx[1][k] + dvx[omp_get_thread_num()][1] * dwdx[0][k];
                    hyy = 2.f * dvx[omp_get_thread_num()][1] * dwdx[1][k] - dvx[omp_get_thread_num()][0] * dwdx[0][k];
                }

                hxx *= 2.f / 3.f;
                hyy *= 2.f / 3.f;

                if (parameters::dim == 1)
                {
                    txx[i] += particles[j].mass * hxx / particles[j].rho;
                    txx[j] += particles[i].mass * hxx / particles[i].rho;
                }
                else if (parameters::dim == 2)
                {
                    txx[i] += particles[j].mass * hxx / particles[j].rho;
                    txx[j] += particles[i].mass * hxx / particles[i].rho;

                    txy[i] += particles[j].mass * hxy / particles[j].rho;
                    txy[j] += particles[i].mass * hxy / particles[i].rho;

                    tyy[i] += particles[j].mass * hyy / particles[j].rho;
                    tyy[j] += particles[i].mass * hyy / particles[i].rho;

                }
            }
        }

#pragma omp for private(pressure) schedule(static)  
        for (int i = 0; i < n_total; i++)
        {
            if (parameters::visc)
            {
                if (parameters::dim == 1)
                {
                    tdsdt[i] = txx[i] * txx[i];
                }
                else if (parameters::dim == 2)
                {
                    tdsdt[i] = txx[i] * txx[i] + 2.f * txy[i] * txy[i] + tyy[i] * tyy[i];
                }
                tdsdt[i] = 0.5f * eta[i] / particles[i].rho * tdsdt[i];
            }

            if (abs(particles[i].type) == 1)  pressure.p_gas(particles[i].rho, particles[i].u, &particles[i].p, &particles[i].c);
            else if (abs(particles[i].type) == 2) pressure.p_art_water(particles[i].rho, &particles[i].p, &particles[i].c);
        }

#pragma omp for private(h, he, rhoij) schedule(static)  
        for (int k = 0; k < niac; k++)
        {
            int i = pair_i[k];
            int j = pair_j[k];
            he = 0.f;

            rhoij = 1.0f / (particles[i].rho * particles[j].rho);
            //  std::cout << rhoij << std::endl;
            if (parameters::paSph == 1)
            {
                for (int d = 0; d < parameters::dim; d++)
                {
                    h = -(particles[i].p + particles[j].p) * dwdx[d][k];
                    he = he + (particles[j].vx[d] - particles[i].vx[d]) * h;
                    if (parameters::visc)
                    {
                        if (d == 0)
                        {
                            h = h + (eta[i] * txx[i] + eta[j] * txx[j]) * dwdx[0][k];
                            if (parameters::dim >= 2)  h = h + (eta[i] * txy[i] + eta[j] * txy[j]) * dwdx[1][k];
                        }
                        else if (d == 1) h = h + (eta[i] * txy[i] + eta[j] * txy[j]) * dwdx[0][k] + (eta[i] * tyy[i] + eta[j] * tyy[j]) * dwdx[1][k];
                    }
                    h = h * rhoij;

                    dvxdt[d][i] = dvxdt[d][i] + particles[j].mass * h;
                    dvxdt[d][j] = dvxdt[d][j] - particles[i].mass * h;
                }
                he = he * rhoij;
                dedt[i] = dedt[i] + particles[j].mass * he;
                dedt[j] = dedt[j] + particles[i].mass * he;
            }
            else if (parameters::paSph == 2)
            {
                for (int d = 0; d < parameters::dim; d++)
                {
                    h = -(particles[i].p / particles[i].rho / particles[i].rho + particles[j].p / particles[j].rho / particles[j].rho) * dwdx[d][k];

                    he += (particles[j].vx[d] - particles[i].vx[d]) * h;
                    //std::cout <<he << std::endl;

                    if (parameters::visc)
                    {
                        if (d == 0)
                        {
                            h += (eta[i] * txx[i] / particles[i].rho / particles[i].rho
                                + eta[j] * txx[j] / particles[j].rho / particles[j].rho) * dwdx[0][k];
                            if (parameters::dim >= 2)
                            {
                                h += (eta[i] * txy[i] / particles[i].rho / particles[i].rho
                                    + eta[j] * txy[j] / particles[j].rho / particles[j].rho) * dwdx[1][k];
                            }
                        }
                        else if (d == 1)
                        {
                            h = h + (eta[i] * txy[i] / particles[i].rho / particles[i].rho
                                + eta[j] * txy[j] / particles[j].rho / particles[j].rho) * dwdx[0][k]
                                + (eta[i] * tyy[i] / particles[i].rho / particles[i].rho
                                    + eta[j] * tyy[j] / particles[j].rho / particles[j].rho) * dwdx[1][k];
                        }
                    }

                    dvxdt[d][i] = dvxdt[d][i] + particles[j].mass * h;
                    dvxdt[d][j] = dvxdt[d][j] - particles[i].mass * h;
                }

                dedt[i] = dedt[i] + particles[j].mass * he;
                dedt[j] = dedt[j] + particles[i].mass * he;
            }
        } 
#pragma omp for  schedule(static)  
        for (int i = 0; i < n_total; i++)
        {
            dedt[i] = tdsdt[i] + 0.5f * dedt[i];
        }
    }
}