#include <iostream>
#include "SPH_Method.h"

int main()
{
	omp_set_num_threads(parameters::num_thread);
	std::string path = "C:\\Users\\Ilya\\Desktop\\Save\\";
	double dt = 0.005;
	int maxstep = 0;
	SPH sph;
	if (parameters::shockTube) 	dt = 0.005;
	if (parameters::shearCavity)  dt = 5.0e-5;
	if (parameters::dambBrake)  dt = 2.8e-6;


	maxstep = 1000;
	sph.loadParticle();

	std::cout << "Enter maximal time steps: ";
	std::cin >> maxstep;
	sph.timeIntegration(dt, maxstep);

	sph.saveParticle(path + "f_xv -", path + "f_state -", path + "f_other -", std::to_string(maxstep));

	return 0;
}
