#include"neighbor_search.h"
SearchNNSP::SearchNNSP()
{
    dxiac_parallel = new double* [parameters::num_thread];
    tdwdx_parallel = new double*[parameters::dim];
    for (int i = 0; i < parameters::num_thread; i++)
    {
        dxiac_parallel[i] = new double[parameters::dim];
        tdwdx_parallel[i] = new double[parameters::dim];
    }
    dxiacFG = new double[parameters::dim];
    tdwdxFG = new double[parameters::dim];
    _niac_parallel = new int[parameters::num_thread];
}

SearchNNSP::~SearchNNSP()
{
    for (int i = 0; i < parameters::num_thread; i++)
    {
        delete[] tdwdx_parallel[i];
        delete[] dxiac_parallel[i];
    }
    delete[] dxiac_parallel;
    delete[] tdwdx_parallel;

    delete[] dxiacFG;
    delete[] tdwdxFG;

    //delete[] dxiac;
    //delete[] tdwdx;
    delete[] _niac_parallel;
}
void SearchNNSP::SimpleSearch(int n_total, Particle* particles, int* niac, int* pair_i, int* pair_j, double* w, double** dwdx)
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
double SearchNNSP::max_hsml(int  n_total,Particle* particles)
{
    double max = particles[0].hsml;
    for (int i = 1; i < n_total; i++)
    {
        if (max < particles[i].hsml)  max = particles[i].hsml;
    }
    return max;
}
void SearchNNSP::Find_Grid(int n_total, Particle* particles, int* niac, int* pair_i, int* pair_j, double* w, double** dwdx)
{
    double r, driac, mhsml,tmp_w_d;
    int _niac = 0;
    int indexX, indexY;
    int sizeY = 1;
    double hsml = particles[0].hsml;

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
//#pragma omp parallel for private(r,driac,tmp_w_d,indexX,indexY,mhsml)
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

                dxiacFG[0] = particles[i].x[0] - particles[j].x[0];
                driac = dxiacFG[0] * dxiacFG[0];

                for (int d = 1; d < parameters::dim; d++)
                {
                    dxiacFG[d] = particles[i].x[d] - particles[j].x[d];
                    driac = driac + dxiacFG[d] * dxiacFG[d];
                }
                mhsml = (particles[i].hsml + particles[j].hsml) / 2.0f;

                if (sqrt(driac) < (NeighbourList::GetScaleK() * mhsml))
                {
                    if (_niac < parameters::maxInteration)
                    {
                        
                        pair_i[_niac] = i;
                        pair_j[_niac] = j;
                        r = sqrt(driac);
                        kernel.findKernel(r, dxiacFG, hsml, &w[_niac], tdwdxFG);
                        for (int d = 0; d < parameters::dim; d++)
                        {
                            dwdx[d][_niac] = tdwdxFG[d];
                        }

                        _niac++;
                    }
                }
            }
        }
    }
    *niac = _niac;
}

void make_grid(int n_total, Particle* particles, double hsml,int indexX,int indexY, NeighbourList::ParticleGrid grid)
{
    for (int i = 0; i < n_total; i++)
    {
        if (particles[i].type == 0) continue;
        indexX = NeighbourList::Block_x(particles[i].x, hsml);
        indexY = 0;

        if (parameters::dim == 2) indexY = NeighbourList::Block_y(particles[i].x, hsml);
        grid.AddParticles(indexX, indexY, i);
    }
}

void SearchNNSP::parallelFindGrid(int n_total, Particle* particles, int* niac, int* pair_i, int* pair_j, double* w, double** dwdx)
{
    double r, driac, mhsml, tmp_w_d;
    int _niac = 0;
    int indexX =0, indexY =0;
    int sizeY = 1;
    double hsml = particles[0].hsml;

    std::vector<std::vector<double>> tmp_pair_i(parameters::num_thread);
    std::vector<std::vector<double>> tmp_pair_j(parameters::num_thread);
    std::vector<std::vector<double>> tmp_w(parameters::num_thread);
    std::vector<std::vector<std::vector<double>>> tmp_tdwdx(parameters::num_thread);

    //std::vector<std::vector<double>> aa(parameters::num_thread,(parameters::dim,std::vector));
    
    tmp_pair_i.clear();
    tmp_pair_j.clear();
    tmp_w.clear();
    tmp_tdwdx.clear();

    if (parameters::dim == 2) sizeY = NeighbourList::SizeY(hsml);
    NeighbourList::ParticleGrid grid(NeighbourList::SizeX(hsml), sizeY);
    grid.Clear();

    make_grid(n_total,particles,hsml,indexX,indexY,grid);
   
    //#pragma omp parallel for private(r,driac,tmp_w_d,indexX,indexY,mhsml)
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

                dxiac_parallel[omp_get_thread_num()][0] = particles[i].x[0] - particles[j].x[0];
                driac = dxiac_parallel[omp_get_thread_num()][0] * dxiac_parallel[omp_get_thread_num()][0];

                for (int d = 1; d < parameters::dim; d++)
                {
                    dxiac_parallel[omp_get_thread_num()][d] = particles[i].x[d] - particles[j].x[d];
                    driac = driac + dxiac_parallel[omp_get_thread_num()][d] * dxiac_parallel[omp_get_thread_num()][d];
                }
                mhsml = (particles[i].hsml + particles[j].hsml) / 2.0f;

                if (sqrt(driac) < (NeighbourList::GetScaleK() * mhsml))
                {
                    if (_niac < parameters::maxInteration)
                    {
                        tmp_pair_i[omp_get_thread_num()].push_back(i);
                        tmp_pair_j[omp_get_thread_num()].push_back(j);
                       /* pair_i[_niac] = i;
                        pair_j[_niac] = j;*/
                        r = sqrt(driac);
                        kernel.findKernel(r, dxiac_parallel[omp_get_thread_num()], hsml, &tmp_w_d, tdwdx_parallel[omp_get_thread_num()]);
                        tmp_w[omp_get_thread_num()].push_back(tmp_w_d);
                        
                        
                        tmp_tdwdx[omp_get_thread_num()].push_back((std::vector<double>()));
                        for (int d = 0; d < parameters::dim; d++)
                        {
                            tmp_tdwdx[omp_get_thread_num()].push_back((std::vector<double>()));

                           // dwdx[d][_niac] = tdwdxFG[d];
                        }

                        _niac++;
                    }
                }
            }
        }
    }
    *niac = _niac;
}