#pragma once
#include <omp.h>
#include <string>
#include <fstream>
namespace parameters
{
    //dimension
    extern  int dim;
    // max particles
    extern int maxn;
    //max iteration
    extern  int maxInteration;
    // number of threads
    extern int num_thread;

    extern  int paSph;//1
    extern   int nnps;//1 - simple, 2 - list, 3 - tree algorithm
    extern  int sle; // 0 - unchanged 1,2,3
    extern   int skf;
    extern   bool summationDensity,
        averageVelocity,
        virtualPart,
        vpInput,
        visc,
        exForce,
        viscArtificial,
        heatArtificial,
        selfGravity,
        norDensity;
    extern   bool confingInput;
    extern  int nSym; // no symmetry/axis symmetry/center symmetry
    extern   bool intStat;
    extern   int printStep,
        saveStep,
        moniParticle;
    extern  bool shockTube, shearCavity, example, dambBrake;

    extern double xMaxGeom;
    extern double  xMinGeomy;
    extern  double yMaxGeom;
    extern double yMinGeomy;
    extern double zMaxGeom;
    extern double zMinGeomy;

    extern double diff;


    extern double etq;
    extern double alpha;
    extern double beta;

    extern double epsilon;
    extern double rr0;
    extern double dd;
    extern double p1;
    extern double p2;


    extern double gamma_water;
    extern double rho0;
    extern double c0;
    extern double gamma_gas;

    extern double visc_water;
    extern double visc_gas;

    extern double q1;
    extern double q2;

    extern double step_time;
    extern int max_step;

    extern std::string path_Param;
    extern std::string path_Coef;
    extern std::string _path;
    extern std::string _pathSave;

    extern bool compress_Morris;


    void set_params( std::string path_coef, std::string path_param, std::string path, std::string pathSave);
}