#include <iostream>
#include "SPH_Method.h"
#include <string>
int main(int argc, char* argv[])
{
	SPH sph;
	//std::string path = "C:\\Users\\Ilya\\Desktop\\Save\\";
	double dt = 0.005;
	int maxstep = 0;
	omp_set_num_threads(parameters::num_thread);
	double t1, time;
   // parameters::_path = "C:\\Users\\Ilya\\Desktop\\SPHExp\\DambBreak\\Input\\";
    //parameters::_pathSave= "C:\\Users\\Ilya\\Desktop\\SPHExp\\DambBreak\\Save\\";
	if (argc > 1) 
	{
#pragma region set_param
        parameters::_path = argv[3];
        parameters::_pathSave = argv[4];
        std::cout << parameters::_pathSave << std::endl;
        std::ifstream sr1(argv[2]);
        std::ifstream sr2(argv[1]);

        if (sr1.is_open())
        {
            sr1 >> parameters::skf;
            sr1 >> parameters::sle;
            sr1 >> parameters::nnps;
            sr1 >> parameters::paSph;
            sr1 >> parameters::averageVelocity;
            sr1 >> parameters::norDensity;
            sr1 >> parameters::visc;
            sr1 >> parameters::viscArtificial;
            sr1 >> parameters::heatArtificial;
            sr1 >> parameters::confingInput;
            sr1 >> parameters::vpInput;
            sr1 >> parameters::virtualPart;
            sr1 >> parameters::exForce;
            sr1 >> parameters::selfGravity;
            sr1 >> parameters::step_time;
            sr1 >> parameters::max_step;
            sr1 >> parameters::saveStep;
            sr1 >> parameters::dim;
            sr1 >> parameters::xMaxGeom;
            sr1 >> parameters::xMinGeomy;
            sr1 >> parameters::yMaxGeom;
            sr1 >> parameters::yMinGeomy;
            sr1 >> parameters::maxn;
            sr1 >> parameters::shockTube;
            sr1 >> parameters::shearCavity;
            sr1 >> parameters::compress_Morris;

            parameters::maxInteration = parameters::maxn * 100;
        }
        else
        {
            std::cout << "No find file" << std::endl;
        }

        if(sr2.is_open()) 
        {
            sr2 >> parameters::etq;
            sr2 >> parameters::alpha;
            sr2 >> parameters::beta;

            sr2 >> parameters::epsilon;

            sr2 >> parameters::rr0;
            sr2 >> parameters::dd;
            sr2 >> parameters::p1;
            sr2 >> parameters::p2;

            sr2 >> parameters::gamma_water;
            sr2 >> parameters::rho0;
            sr2 >> parameters::c0;
            sr2 >> parameters::gamma_gas;

            sr2 >> parameters::visc_water;
            sr2 >> parameters::visc_gas;

            sr2 >> parameters::q1;
            sr2 >> parameters::q2;
        }

       


        sr1.close();
        sr2.close();

        maxstep = parameters::max_step;
        dt = parameters::step_time;
#pragma endregion

        std::cout << "In this simulation, the time step: " << dt << std::endl
            << "and the maximum number of steps are used: " << maxstep << std::endl;
	}
	else
	{
		
		std::cout << "Enter maximal time steps: ";
		std::cin >> maxstep;
		
		if (parameters::shockTube) 	dt = 0.005;
		if (parameters::shearCavity)  dt = 5.0e-5;
		if (parameters::dambBrake)  dt = 5.e-6;
	}

	sph.loadParticle();

	
	t1 = omp_get_wtime();
	sph.timeIntegration(dt, maxstep);
	sph.saveParticle(parameters::_pathSave + "f_xv -", parameters::_pathSave + "f_state -", parameters::_pathSave + "f_other -", std::to_string(maxstep));
	std::cout << omp_get_wtime() - t1 << std::endl;
	system("pause");
	return 0;
}
