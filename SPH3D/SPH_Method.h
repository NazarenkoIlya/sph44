#pragma once
#include "particles.h"
#include "param.h"
#include "linked_list.h"
#include "art_heat.h"  
#include "art_visc.h"
#include "density.h"
#include "ext_force.h"
#include "kernel.h"
#include "int_force.h"
#include "hsml.h"
#include "input.h"
#include "input_virtual.h"
#include "output.h"
//#include <iostream>
#include <fstream>
#include <string>
#include<array>
class SPH
{

private: 
    int itimestep = 0, nstart = 0, n_total = 0;
    double time = 0;
    double** v_min, ** dx, * rho_min, * u_min, ** av, ** dvx, * du, * tdsdt, * drho;
    Particle* particles;
    ArtHeat art_heat;
    ArtVisc art_visc;
    Density density;
    ExtForce ext_force;
    HSML hsml;
    IntForce int_force;
    Kernel kernel;
    InitVirtParticle initVirtParticle;
   
    SavingParticle saveParticles;


    double* findViscosity(int ntotal, double* eta);
    void stepSingle(double dt, int nTotal, Particle* particles, int itimestep, double* tdsdt, double** dx, double** dvx, double* du, double** av, double* drho);
    double** findAvVel(Particle* particles, int niac, int* pair_i, int* pair_j, double* w, double** av);
    void SimpleSearch(int n_total, Particle* particles, int* niac, int* pair_i, int* pair_j, double* w, double** dwdx);
    double max_hsml();
    void Find_Grid(int n_total, int* niac, int* pair_i, int* pair_j, double* w, double** dwdx);

public:
    SPH();
    void timeIntegration(double dt, int maxtimestep);
    void loadParticle();
    void saveParticle(std::string name_file1, std::string name_file2, std::string name_file3, std::string iter);
    ~SPH();
};
