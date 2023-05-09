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
#include "av_vel.h"
#include "output.h"
#include "neighbor_search.h"
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
    SearchNNSP nnsp;
    //Kernel kernel;
    InitVirtParticle initVirtParticle;
    SavingParticle saveParticles;
    AvVel av_vel;

    double* findViscosity(int ntotal, double* eta);
    void stepSingle(double dt, int nTotal, Particle* particles, int itimestep, double* tdsdt, double** dx, double** dvx, double* du, double** av, double* drho);
   
public:
    SPH();
    void timeIntegration(double dt, int maxtimestep);
    void loadParticle();
    void saveParticle(std::string name_file1, std::string name_file2, std::string name_file3, std::string iter);
    ~SPH();
};
