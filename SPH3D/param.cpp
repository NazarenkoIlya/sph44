#include "param.h"
int parameters::dim = 2;
int parameters::maxn = 100000;
int parameters::maxInteration = maxn * 100;
int parameters::num_thread = 1;


double parameters::xMaxGeom = 0.003e0;
double parameters::xMinGeomy = 0.e0;
double parameters::yMaxGeom = 0.003e0;
double parameters::yMinGeomy = 0.e0;

//double parameters::xMaxGeom = 5.e0;
//double  parameters::xMinGeomy = -5.e0;
//double parameters::yMaxGeom = 2.e0;
//double parameters::yMinGeomy = -2.e0;

double parameters::zMaxGeom = 2.e0;
double parameters::zMinGeomy = -2.e0;

int parameters::paSph = 2;//1
int parameters::nnps = 2;//1 - simple, 2 - list, 3 - tree algorithm
int parameters::sle = 0; // 0 - unchanged 1,2,3
int parameters::skf = 1;

bool parameters::summationDensity = true;
bool parameters::averageVelocity = false;
bool parameters::confingInput = false;
bool parameters::virtualPart = false;
bool parameters::vpInput = false;
bool parameters::visc = false;
bool parameters::exForce = false;
bool parameters::viscArtificial = true;
bool parameters::heatArtificial = false;
bool parameters::selfGravity = false;
bool parameters::norDensity = false;

int parameters::nSym = 0; // no symmetry/axis symmetry/center symmetry
bool parameters::intStat = true;
int parameters::printStep = 100;
int parameters::saveStep = 200;
int parameters::moniParticle = 1600;

bool parameters::shockTube = false;
bool parameters::shearCavity = false;
bool parameters::dambBrake = true;
bool parameters::example = false;

double parameters::diff = 0.;

double parameters::etq = 0.1;
double parameters::alpha = 1.;
double parameters::beta = 1.;

double parameters::epsilon = 0.3;

double parameters::rr0 = 1.25e-5;
double parameters::dd = 1.e-2;
double parameters::p1 = 12;
double parameters::p2 = 4;

double parameters::gamma_water = 7;
double parameters::rho0 = 1000;
double parameters::c0 = 0.1;
double parameters::gamma_gas = 1.4;

double parameters::visc_water = 1.e-3;
double parameters::visc_gas = 0;

double parameters::q1 = 0.1;
double parameters::q2 = 1.;

double parameters::step_time = 0;
double parameters::max_step = 0;
std::string parameters::_path = "C:\\Users\\Ilya\\Desktop\\Input\\";
std::string parameters::_pathSave = "C:\\Users\\Ilya\\Desktop\\Save\\";

void parameters::set_params(std::string path_param, std::string path_coef,std::string path, std::string pathSave)
{
    _path = path;
    _pathSave = pathSave;
    std::ifstream sr1(path_param);
    std::ifstream sr2(path_coef);

    sr1 >> skf;
    sr1 >> sle;
    sr1 >> nnps;
    sr1 >> paSph;
    sr1 >> averageVelocity;
    sr1 >> norDensity;
    sr1 >> visc;
    sr1 >> viscArtificial;
    sr1 >> heatArtificial;
    sr1 >> confingInput;
    sr1 >> vpInput;
    sr1 >> virtualPart;
    sr1 >> exForce;
    sr1 >> selfGravity;
    sr1 >> step_time;
    sr1 >> max_step;
    sr1 >> saveStep;
    sr1 >> dim;
    sr1 >> xMaxGeom;
    sr1 >> xMinGeomy;
    sr1 >> yMaxGeom;
    sr1 >> yMinGeomy;
    sr1 >> maxn;
    sr1 >> shockTube;
    sr1 >> shearCavity;

    maxInteration = maxn * 100;

    sr2 >> etq;
    sr2 >> alpha;
    sr2 >> beta;

    sr2 >> epsilon;

    sr2 >> rr0;
    sr2 >> dd;
    sr2 >> p1;
    sr2 >> p2;

    sr2 >> gamma_water;
    sr2 >> rho0;
    sr2 >> c0;
    sr2 >> gamma_gas;

    sr2 >> visc_water;
    sr2 >> visc_gas;

    sr2 >> q1;
    sr2 >> q2;


    sr1.close();
    sr2.close();
}
