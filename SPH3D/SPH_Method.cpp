#include"SPH_Method.h"
#include <iostream>
//#include <iostream>

   
    double* SPH::findViscosity(int ntotal, double* eta)
    {
#pragma omp parallel for 
        for (int i = 0; i < ntotal; i++)
        {
            if (abs(particles[i].type) == 1)    eta[i] = parameters::visc_gas;
            else if (abs(particles[i].type) == 2) eta[i] = parameters::visc_water;
        }
        return eta;
    }
    void SPH::stepSingle(double dt, int nTotal, Particle* particles, int itimestep, double* tdsdt, double** dx, double** dvx, double* du, double** av, double* drho)
    {
        double t1, t2,time, tt1,tt2;
        int nVirt, niac = 0;
        int* pair_i = new int[parameters::maxInteration];
        int* pair_j = new int[parameters::maxInteration];
        double* w = new double[parameters::maxInteration];
        int* ns = new int[parameters::maxn];
        double* eta = new double[parameters::maxn];
        double* avdudt = new double[parameters::maxn];
        double* ahdudt = new double[parameters::maxn];
        double* c = new double[parameters::maxn];

        double** dwdx = new double* [parameters::dim];
        double** indvxdt = new double* [parameters::dim];
        double** exdvxdt = new double* [parameters::dim];
        double** ardvxdt = new double* [parameters::dim];
        for (int d = 0; d < parameters::dim; d++)
        {
            dwdx[d] = new double[parameters::maxInteration];
            indvxdt[d] = new double[parameters::maxn];
            exdvxdt[d] = new double[parameters::maxn];
            ardvxdt[d] = new double[parameters::maxn];
        }

        for (int i = 0; i < n_total; i++)
        {
            avdudt[i] = 0.;
            ahdudt[i] = 0.;
            for (int d = 0; d < parameters::dim; d++)
            {
                indvxdt[d][i] = 0.;
                ardvxdt[d][i] = 0.;
                exdvxdt[d][i] = 0.;
            }
        }

        nVirt = 0;

        if (parameters::virtualPart)  initVirtParticle.createVirtPart(itimestep, &nVirt, n_total,particles,"C:\\Users\\Ilya\\Desktop\\Input\\xv_vp.txt", "C:\\Users\\Ilya\\Desktop\\Input\\state_vp.txt", "C:\\Users\\Ilya\\Desktop\\Input\\other_vp.txt");
       
        t1 = omp_get_wtime();
        tt1 = omp_get_wtime();
        if (parameters::nnps == 1)nnsp.SimpleSearch(n_total + nVirt, particles,&niac, pair_i, pair_j, w, dwdx);
        else  nnsp.Find_Grid(n_total + nVirt, particles,&niac, pair_i, pair_j, w, dwdx);
        //nnsp.parallelFindGrid(n_total + nVirt, particles, &niac, pair_i, pair_j, w, dwdx);
        t2 = omp_get_wtime();
        std::cout << t2 - t1 << std::endl;
        
        t1 = omp_get_wtime();
        if (parameters::summationDensity) density.findSumDensity(n_total + nVirt, niac, particles, pair_i, pair_j, w);
        else drho = density.findConDensity(n_total + nVirt, particles, niac, pair_i, pair_j, dwdx, drho);
        t2 = omp_get_wtime();
        std::cout << t2 - t1 << std::endl;
        
        t1 = omp_get_wtime();
        if (parameters::visc) eta = findViscosity(n_total + nVirt, eta);
        t2 = omp_get_wtime();
        std::cout << t2 - t1 << std::endl;
       
        t1 = omp_get_wtime();
        int_force.find(n_total + nVirt, particles, niac, pair_i, pair_j, dwdx, indvxdt, tdsdt, du, eta);
        t2 = omp_get_wtime();
        std::cout << t2 - t1 << std::endl;

        t1 = omp_get_wtime();
        if (parameters::viscArtificial)  art_visc.find(n_total + nVirt, particles, niac, pair_i, pair_j, dwdx, ardvxdt, avdudt);
        t2 = omp_get_wtime();
        std::cout << t2 - t1 << std::endl;

        t1 = omp_get_wtime();
        if (parameters::exForce) exdvxdt = ext_force.findExtForce(n_total + nVirt, particles, niac, pair_i, pair_j, exdvxdt);
        t2 = omp_get_wtime();
        std::cout<< t2 - t1 << std::endl;

        t1 = omp_get_wtime();
        if (parameters::sle != 0)  hsml.upgradeHsml(dt, n_total, particles, niac, pair_i, pair_j, dwdx);
        t2 = omp_get_wtime();
        std::cout << t2 - t1 << std::endl;
        
        t1 = omp_get_wtime();
        if (parameters::heatArtificial) ahdudt = art_heat.find(n_total + nVirt, niac, pair_i, pair_j, w, dwdx, particles);
        t2 = omp_get_wtime();
        std::cout<< t2 - t1 << std::endl;
        t1 = omp_get_wtime();
        if (parameters::averageVelocity) av = av_vel.findAvVel(particles,n_total,niac, pair_i, pair_j, w, av);
        t2 = omp_get_wtime();
        std::cout << t2 - t1 << std::endl;
        std::cout << std::endl;
#pragma  omp parallel for schedule(static)
        for (int i = 0; i < n_total; i++)
        {
            for (int d = 0; d < parameters::dim; d++)
            {
                dvx[d][i] = indvxdt[d][i] + exdvxdt[d][i] + ardvxdt[d][i];
            }
            du[i] = du[i] + avdudt[i] + ahdudt[i];
        }  
        
        tt2 = omp_get_wtime();
        std::cout << tt2 - tt1 << std::endl;
        std::cout << std::endl;
        delete[] pair_i;
        delete[] pair_j;
        delete[] ns;
        delete[] w;
        delete[] avdudt;
        delete[] ahdudt;
        delete[] c;

        for (int d = 0; d < parameters::dim; d++)
        {
            delete[] dwdx[d];
        }
        delete[] dwdx;

        for (int d = 0; d < parameters::dim; d++)
        {
            delete[] indvxdt[d];
        }
        delete[] indvxdt;

        for (int d = 0; d < parameters::dim; d++)
        {
            delete[] exdvxdt[d];
        }
        delete[] exdvxdt;

        for (int d = 0; d < parameters::dim; d++)
        {
            delete[] ardvxdt[d];
        }
        delete[] ardvxdt;
        delete[] eta;
    }
    
   
    
   
    SPH::SPH()
    {
       
        du = new double[parameters::maxn];
        u_min = new double[parameters::maxn];
        rho_min = new double[parameters::maxn];
        tdsdt = new double[parameters::maxn];
        drho = new double[parameters::maxn];

        v_min = new double* [parameters::dim];
        dx = new double* [parameters::dim];
        av = new double* [parameters::dim];
        dvx = new double* [parameters::dim];
        for (int d = 0; d < parameters::dim; d++)
        {
            v_min[d] = new double[parameters::maxn];
            dx[d] = new double[parameters::maxn];
            av[d] = new double[parameters::maxn];
            dvx[d] = new double[parameters::maxn];
        }

        for (int i = 0; i < parameters::maxn; i++)
        {
            u_min[i] = 0.f;
            du[i] = 0.f;
            rho_min[i] = 0.f;
            tdsdt[i] = 0.f;
            for (int d = 0; d < parameters::dim; d++)
            {
                v_min[d][i] = 0.;
                av[d][i] = 0.;
                dvx[d][i] = 0.;
                dx[d][i] = 0.;
            }
        }
        particles = new Particle[parameters::maxn];
        
      
    }
    void SPH::timeIntegration(double dt, int maxtimestep)
    {
        for (int i = 0; i < n_total; i++)
        {
            for (int d = 0; d < parameters::dim; d++)
            {
                av[d][i] = 0.f;
            }
        }
        double temp_u = 0, temp_rho = 0;

        for (itimestep = 0; itimestep < maxtimestep; itimestep++)
        {        
            if (itimestep != 0) 
            {
#pragma  omp parallel  for private(temp_u,temp_rho) 
                for (int i = 0; i < n_total; i++)
                {
                    u_min[i] = particles[i].u;
                    temp_u = 0.;

                    if (parameters::dim == 1) temp_u = -parameters::nSym * particles[i].p * particles[i].vx[0] / particles[i].x[0] / particles[i].rho;
                   
                    particles[i].u += ((dt / 2.0) * (du[i] + temp_u));
                    
                    if (particles[i].u < 0) particles[i].u = 0.;

                    if (!parameters::summationDensity)
                    {
                        rho_min[i] = particles[i].rho;
                        temp_rho = 0.;
                        if (parameters::dim == 1)   temp_rho = -parameters::nSym * particles[i].rho * particles[i].vx[0] / particles[i].x[0];
                        particles[i].rho += (dt / 2.) * (drho[i] + temp_rho);
                    }

                    for (int d = 0; d < parameters::dim; d++)
                    {
                        v_min[d][i] = particles[i].vx[d];
                     
                        particles[i].vx[d] += (dt / 2.) * dvx[d][i];
                       // std::cout << i << "----" << v_min[d][i] << " " << particles[i].vx[d] << std::endl;
                    }
                }
            }
            stepSingle(dt, n_total,particles, itimestep, tdsdt,dx,dvx,du,av,drho);
          
            if (itimestep == 0) 
            {
               
#pragma  omp parallel for private(temp_u,temp_rho)  
                for (int i = 0; i < n_total; i++)
                {
                    temp_u = 0.;
                    if (parameters::dim == 1) temp_u = -parameters::nSym * particles[i].p * particles[i].vx[0] / particles[i].x[0] / particles[i].rho;
                    particles[i].u += ((dt / 2.0) * (du[i] + temp_u));

                    if (particles[i].u < 0)  particles[i].u = 0;

                    if (!parameters::summationDensity)
                    {
                        temp_rho = 0.;
                        if (parameters::dim == 1) temp_rho = -parameters::nSym * particles[i].rho * particles[i].vx[0] / particles[i].x[0];
                        particles[i].rho += ((dt / 2.0) * (du[i] + temp_u));
                    }

                    for (int d = 0; d < parameters::dim; d++)
                    {
                        particles[i].vx[d] += ((dt / 2.) * dvx[d][i] + av[d][i]);
                        particles[i].x[d] += dt * particles[i].vx[d];
                    }
                }

            }
            else 
            {
#pragma  omp parallel for private(temp_u,temp_rho)
                for (int i = 0; i < n_total; i++)
                {
                    temp_u = 0.;
                    if (parameters::dim == 1) temp_u = -parameters::nSym * particles[i].p * particles[i].vx[0] / particles[i].x[0] / particles[i].rho;

                    particles[i].u = u_min[i] + dt * (du[i] + temp_u);
                    if (particles[i].u < 0) particles[i].u = 0.;

                    if (!parameters::summationDensity)
                    {
                        temp_rho = 0.;
                        if (parameters::dim == 1) temp_rho = -parameters::nSym * particles[i].rho * particles[i].vx[0] / particles[i].x[0];
                        particles[i].rho = rho_min[i] + dt * (drho[i] + temp_rho);
                    }
                    for (int d = 0; d < parameters::dim; d++)
                    {
                        particles[i].vx[d] = v_min[d][i] + dt * dvx[d][i] + av[d][i];
                        particles[i].x[d] += dt * particles[i].vx[d];
                    }
                }
            }  

            time += dt;
            if ((itimestep % parameters::saveStep) == 0)
            {
                std::string path = "C:\\Users\\Ilya\\Desktop\\Save\\";
                saveParticle(path + "f_xv -", path + "f_state -", path + "f_other -",std::to_string(itimestep));
            }
            
        }
    }

    void SPH::loadParticle()
    {
        ParticleInitialization particlesInit;
        particlesInit.loadInitialParticle(particles, &n_total);
    }
    void SPH::saveParticle(std::string name_file1, std::string name_file2, std::string name_file3, std::string iter)
    {
        std::string path = "C:\\Users\\Ilya\\Desktop\\Save\\";
        saveParticles.save(n_total,particles,path + "f_xv -", path + "f_state -", path + "f_other -", std::to_string(itimestep));
    }
    SPH::~SPH()
    {
        for (int d = 0; d < parameters::dim; d++)
        {
            delete[] v_min[d];
        }
        delete[] v_min;

        for (int d = 0; d < parameters::dim; d++)
        {
            delete[] av[d];
        }
        delete[] av;

        for (int d = 0; d < parameters::dim; d++)
        {
            delete[] dvx[d];
        }
        delete[] dvx;

        for (int d = 0; d < parameters::dim; d++)
        {
            delete[] dx[d];
        }
        delete[] dx;

        delete[] du;
        delete[] u_min;
        delete[] rho_min;
        delete[] tdsdt;
        delete[] drho;
        delete[] particles;
    }
