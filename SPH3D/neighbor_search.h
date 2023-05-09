#pragma once
#include "param.h"
#include "particles.h"
#include "kernel.h"
#include "linked_list.h"
class SearchNNSP
{
public:
	void SimpleSearch(int n_total, Particle* particles, int* niac, int* pair_i, int* pair_j, double* w, double** dwdx);
	void Find_Grid(int n_total, Particle* particles,int* niac, int* pair_i, int* pair_j, double* w, double** dwdx);
	void parallelFindGrid(int n_total, Particle* particles, int* niac, int* pair_i, int* pair_j, double* w, double** dwdx);
	SearchNNSP();
	~SearchNNSP();

private:
	double max_hsml(int n_total, Particle* particles);
	Kernel kernel;
	double* dxiacFG;
	double* tdwdxFG;

	double** dxiac_parallel;
	double** tdwdx_parallel;

	int* _niac_parallel;
};

