#include "ext_force.h"
#include <iostream>

ExtForce::ExtForce()
{
    dx = new double[parameters::dim];
    //rr0 = 4.9e-5f;
    //dd = 1.e-5f;
  
    rr0 = 1.25e-5f;
    dd = 1.e-2f;
    p1 = 12.f;
    p2 = 6.f;

}

ExtForce::~ExtForce()
{
    delete[] dx;
}
//double max_velocity(Particle* particles, int n_total)
//{
//    double max = -10;
//    for (int i = 0; i < n_total; i++)
//    {
//        for (int d = 0; d < parameters::dim; d++)
//        {
//            double sqr = particles[i].vx[d] * particles[i].vx[d];
//            if (max < sqr)
//            {
//                max = sqr;
//            }
//        }
//    }
//    return max;
//}

 double** ExtForce::findExtForce(int n_total, Particle* particles, int niac, int* pair_i, int* pair_j, double** dvxdt)
{

     //if (parameters::dambBrake)
     //{
     //    rr0 = parameters::diff;
     //    dd = max_velocity(particles, n_total);
     //}
     double rr, f;
//#pragma omp parallel for
    for (int i = 0; i < n_total; i++)
    {
        for (int d = 0; d < parameters::dim; d++)
        {
            dvxdt[d][i] = 0.f;
        }
    }

    if (parameters::selfGravity)
    {
//#pragma omp parallel for
        for (int i = 0; i < n_total; i++)
        {
            dvxdt[parameters::dim - 1][i] = -9.8f;
        }
    }

//#pragma omp parallel for private (rr,f)
    for (int k = 0; k < niac; k++)
    {
        int i = pair_i[k];
        int j = pair_j[k];

        if (particles[i].type > 0 && particles[j].type < 0)
        {
            rr = 0.;
            for (int d = 0; d < parameters::dim; d++)
            {
                dx[d] = particles[i].x[d] - particles[j].x[d];
                rr += dx[d] * dx[d];
            }
            rr = sqrt(rr);
            if (rr < rr0)
            {
                f = (pow(rr0 / rr, p1) - pow(rr0 / rr, p2)) / (rr * rr);
                for (int d = 0; d < parameters::dim; d++)
                {
                    dvxdt[d][i] += dx[d] * f * dd;
                }
            }
        }
    }
    return dvxdt;
}