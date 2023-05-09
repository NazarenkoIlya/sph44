#include <iostream>
#include "SPH_Method.h"

int main(int argc, char* argv[])
{
	omp_set_num_threads(parameters::num_thread);
	double t1, time;
	//std::string path = "C:\\Users\\Ilya\\Desktop\\Save\\";
	double dt = 0.005;
	int maxstep = 0;
	SPH sph;
	if (parameters::shockTube) 	dt = 0.005;
	if (parameters::shearCavity)  dt = 5.0e-5;
	if (parameters::dambBrake)  dt = 5.e-6;


	maxstep = 1000;
	sph.loadParticle();

	std::cout << "Enter maximal time steps: ";
	std::cin >> maxstep;
	t1 = omp_get_wtime();
	sph.timeIntegration(dt, maxstep);

	sph.saveParticle(parameters::_pathSave + "f_xv -", parameters::_pathSave + "f_state -", parameters::_pathSave + "f_other -", std::to_string(maxstep));
	std::cout << omp_get_wtime() - t1 << std::endl;
	return 0;
}
