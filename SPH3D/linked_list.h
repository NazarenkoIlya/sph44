#pragma once
#include "param.h"
//#include <iostream>
#include <vector>
#include <cmath>

namespace NeighbourList
{
	class Particles
	{
	public:
		void operator=(int index_particles);
		std::vector<int> GetParticle();
		void Clear();
	private:
		std::vector<int> index_particle;
	};

	class ParticleGrid
	{
	public:
		ParticleGrid(int n_cells_x, int n_cells_y);
		void AddParticles(int index_cell_x, int index_cell_y, int iParticles);
		void Clear();
		std::vector<Particles> GetNeighbourGridParticles(int index_cell_x, int index_cell_y);
		~ParticleGrid();

	private:
		int n_cells_x;
		int n_cells_y;
		Particles** grid;
		void IndexParticlesOnBoundary(int* index_cell, int size);
		bool isIndexValid2D(int index_cell_x, int index_cell_y);
		bool isIndexValid1D(int index_cell_x);
	};

	double GetScaleK();

	double BlockSize(double hsml);

	int SizeX(double hsml);

	int SizeY(double hsml);

	int Block_x(double* x, double hsml);

	int Block_y(double* x, double hsml);
}