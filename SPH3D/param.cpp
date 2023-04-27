#include "param.h"
int parameters::dim = 2;
int parameters::maxn = 100000;
int parameters::maxInteration = maxn * 100;
int parameters::num_thread = 4;


double parameters::xMaxGeom = 0.002e0;
double parameters::xMinGeomy = 0.e0;
double parameters::yMaxGeom = 0.002e0;
double parameters::yMinGeomy = 0.e0;

//double parameters::xMaxGeom = 5.e0;
//double  parameters::xMinGeomy = -5.e0;
//double parameters::yMaxGeom = 2.e0;
//double parameters::yMinGeomy = -2.e0;

double parameters::zMaxGeom = 2.e0;
double parameters::zMinGeomy = -2.e0;

int parameters::paSph = 2;//1
int parameters::nnps = 1;//1 - simple, 2 - list, 3 - tree algorithm
int parameters::sle = 0; // 0 - unchanged 1,2,3
int parameters::skf = 2;

bool parameters::summationDensity = true;
bool parameters::averageVelocity = true;
bool parameters::confingInput = false;
bool parameters::virtualPart = true;
bool parameters::vpInput = false;
bool parameters::visc = false;
bool parameters::exForce = true;
bool parameters::viscArtificial = false;
bool parameters::heatArtificial = false;
bool parameters::selfGravity = true;
bool parameters::norDensity = true;

int parameters::nSym = 0; // no symmetry/axis symmetry/center symmetry
bool parameters::intStat = true;
int parameters::printStep = 100;
int parameters::saveStep = 1000;
int parameters::moniParticle = 1600;

bool parameters::shockTube = false;
bool parameters::shearCavity = false;
bool parameters::dambBrake = true;
bool parameters::example = false;

double parameters::diff = 0.;