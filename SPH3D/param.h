#pragma once
#include <omp.h>
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
}