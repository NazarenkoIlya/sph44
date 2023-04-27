#include"SPH_Method.h"
#include <iostream>
//#include <iostream>

   
    double* SPH::findViscosity(int ntotal, double* eta)
    {
//#pragma omp parallel for 
        for (int i = 0; i < ntotal; i++)
        {
            if (abs(particles[i].type) == 1)    eta[i] = 0.;
            else if (abs(particles[i].type) == 2) eta[i] = 1.e-3;
        }
        return eta;
    }
    void SPH::stepSingle(double dt, int nTotal, Particle* particles, int itimestep, double* tdsdt, double** dx, double** dvx, double* du, double** av, double* drho)
    {
      
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
       
        if (parameters::nnps == 1) SimpleSearch(n_total + nVirt, particles,&niac, pair_i, pair_j, w, dwdx);
        else  Find_Grid(n_total + nVirt, &niac, pair_i, pair_j, w, dwdx);
        
        if (parameters::summationDensity) density.findSumDensity(n_total + nVirt, niac, particles, pair_i, pair_j, w);
        else drho = density.findConDensity(n_total + nVirt, particles, niac, pair_i, pair_j, dwdx, drho);
   
        if (parameters::visc) eta = findViscosity(n_total + nVirt, eta);
       
        int_force.find(n_total + nVirt, particles, niac, pair_i, pair_j, dwdx, indvxdt, tdsdt, du, eta);
      
        if (parameters::viscArtificial)  art_visc.find(n_total + nVirt, particles, niac, pair_i, pair_j, dwdx, ardvxdt, avdudt);
     
        if (parameters::exForce) exdvxdt = ext_force.findExtForce(n_total + nVirt, particles, niac, pair_i, pair_j, exdvxdt);
       
        if (parameters::sle != 0)  hsml.upgradeHsml(dt, n_total, particles, niac, pair_i, pair_j, dwdx);
       
        if (parameters::heatArtificial) ahdudt = art_heat.find(n_total + nVirt, niac, pair_i, pair_j, w, dwdx, particles);
    
        if (parameters::averageVelocity) av = findAvVel(particles,niac, pair_i, pair_j, w, av);
      #pragma  omp for 
        for (int i = 0; i < n_total; i++)
        {
            for (int d = 0; d < parameters::dim; d++)
            {
                dvx[d][i] = indvxdt[d][i] + exdvxdt[d][i] + ardvxdt[d][i];
              /*  if (i%10==0)
                {
                    std::cout << i << "----\t" << indvxdt[0][i]
                        << "\t" << exdvxdt[0][i]
                        << "\t" << ardvxdt[0][i] << std::endl;
                }
              */

            }
            du[i] = du[i] + avdudt[i] + ahdudt[i];
        }

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
    double** SPH::findAvVel(Particle * particles,int niac, int* pair_i, int* pair_j, double* w, double** av)
    {
        double epsilon = 0.3;
        double* dvx = new double[parameters::dim];;

        for (int i = 0; i < n_total; i++)
        {
            for (int  d = 0; d < parameters::dim; d++)
            {
                av[d][i] = 0.;
            }
        }

//#pragma omp parallel for 
        for (int k = 0; k < niac; k++)
        {
            int i = pair_i[k];
            int j = pair_j[k];

            for (int d = 0; d < parameters::dim; d++)
            {
                dvx[d] = particles[i].vx[d] - particles[j].vx[d];
                av[d][i] = av[d][i] - 2 * particles[j].mass * dvx[d] / (particles[i].rho + particles[j].rho) * w[k];
                av[d][j] = av[d][j] + 2 * particles[i].mass * dvx[d] / (particles[i].rho + particles[j].rho) * w[k];
            }
        }
//#pragma omp parallel for
        for (int i = 0; i < n_total; i++)
        {
            for (int d = 0; d < parameters::dim; d++)
            {
                av[d][i] = epsilon * av[d][i];
            }
        }
        delete[] dvx;
        return av;
    }
    void SPH::SimpleSearch(int n_total, Particle * particles ,int* niac, int* pair_i, int* pair_j, double* w, double** dwdx)
    {
        int scale_k = 2;


        double* dxiac = new double[parameters::dim];
        double* tdwdx = new double[parameters::dim];
        double r, driac, mhsml;
        int i, j, d;
        int _niac = 0;

        if (parameters::skf == 1)  scale_k = 2;
        else if (parameters::skf == 2 || parameters::skf == 3) scale_k = 3;

        for (i = 0; i < n_total - 1; i++)
        {
            for (j = i + 1; j < n_total; j++)
            {
                dxiac[0] = particles[i].x[0] - particles[j].x[0];
                driac = dxiac[0] * dxiac[0];

                for (d = 1; d < parameters::dim; d++)
                {
                    dxiac[d] = particles[i].x[d] - particles[j].x[d];
                    driac = driac + dxiac[d] * dxiac[d];
                }

                mhsml = (particles[i].hsml + particles[j].hsml) / 2.;

                if (sqrt(driac) < (scale_k * mhsml))
                {
                    if (_niac < parameters::maxInteration)
                    {
                        pair_i[_niac] = i;
                        pair_j[_niac] = j;
                        r = sqrt(driac);
                        kernel.findKernel(r, dxiac, mhsml, &w[_niac], tdwdx);
                        for (d = 0; d < parameters::dim; d++)
                        {
                            dwdx[d][_niac] = tdwdx[d];
                        }
                        _niac++;
                    }
                }
            }
        }
        *niac = _niac;
        delete[] dxiac;
        delete[] tdwdx;
    }
    double SPH::max_hsml()
    {
        double max = particles[0].hsml;
        for (int i = 1; i < n_total; i++)
        {
            if (max < particles[i].hsml)  max = particles[i].hsml;
        }
        return max;
    }
    void SPH::Find_Grid(int n_total, int* niac, int* pair_i, int* pair_j, double* w, double** dwdx)
    {
        //std::ofstream sw1("Parallel.txt");
        double* dxiac = new double[parameters::dim];
        double* tdwdx = new double[parameters::dim];
        double r, driac, mhsml;
        int _niac = 0;
        int help = 0;
        int indexX, indexY;
        int sizeY = 1;
        double hsml = max_hsml();

        if (parameters::dim == 2) sizeY = NeighbourList::SizeY(hsml);
        NeighbourList::ParticleGrid grid(NeighbourList::SizeX(hsml), sizeY);
        grid.Clear();

        for (int i = 0; i < n_total; i++)
        {
            if (particles[i].type == 0) continue;

            indexX = NeighbourList::Block_x(particles[i].x, hsml);
            indexY = 0;

            if (parameters::dim == 2) indexY = NeighbourList::Block_y(particles[i].x, hsml);
            grid.AddParticles(indexX, indexY, i);
        }

        for (int i = 0; i < n_total; i++)
        {
            if (particles[i].type == 0) continue;

            indexX = NeighbourList::Block_x(particles[i].x, hsml);
            indexY = 0;

            if (parameters::dim == 2) indexY = NeighbourList::Block_y(particles[i].x, hsml);
            std::vector<NeighbourList::Particles> neighbour_cells;

            neighbour_cells = grid.GetNeighbourGridParticles(indexX, indexY);

            for (auto block : neighbour_cells)
            {
                for (auto j : block.GetParticle())
                {
                    if (i >= j) continue;

                    dxiac[0] = particles[i].x[0] - particles[j].x[0];
                    driac = dxiac[0] * dxiac[0];

                    for (int d = 1; d < parameters::dim; d++)
                    {
                        dxiac[d] = particles[i].x[d] - particles[j].x[d];
                        driac = driac + dxiac[d] * dxiac[d];
                    }
                    mhsml = (particles[i].hsml + particles[j].hsml) / 2.0f;

                    if (sqrt(driac) < (NeighbourList::GetScaleK() * mhsml))
                    {
                        if (_niac < parameters::maxInteration)
                        {
                            pair_i[_niac] = i;
                            pair_j[_niac] = j;
                            r = sqrt(driac);
                            kernel.findKernel(r, dxiac, hsml, &w[_niac], tdwdx);
                            for (int d = 0; d < parameters::dim; d++)
                            {
                                dwdx[d][_niac] = tdwdx[d];
                            }

                            _niac++;
                        }
                    }
                }
            }
        }

        *niac = _niac;

        delete[] dxiac;
        delete[] tdwdx;
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
//#pragma  omp parallel  for private(temp_u,temp_rho) 
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
               
//#pragma  omp parallel for private(temp_u,temp_rho)  
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
//#pragma  omp parallel for private(temp_u,temp_rho)
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
