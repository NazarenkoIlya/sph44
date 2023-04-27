#include "linked_list.h"
//#include <iostream>

void NeighbourList::Particles::operator=(int index_particles)
{
	index_particle.push_back(index_particles);
}


std::vector<int> NeighbourList::Particles::GetParticle()
{
	return index_particle;
}

void NeighbourList::Particles::Clear()
{
	index_particle.clear();
}

NeighbourList::ParticleGrid::ParticleGrid(int n_cells_x, int n_cells_y) :n_cells_x{ n_cells_x }, n_cells_y{ n_cells_y }
{
	grid = new Particles * [n_cells_x];
	for (int i = 0; i < n_cells_x; i++)
	{
		grid[i] = new Particles[n_cells_y];
	}
}


void NeighbourList::ParticleGrid::AddParticles(int index_cell_x, int index_cell_y, int iParticles)
{
	IndexParticlesOnBoundary(&index_cell_x, n_cells_x);
	IndexParticlesOnBoundary(&index_cell_y, n_cells_y);
	grid[index_cell_x][index_cell_y] = iParticles;
}

void NeighbourList::ParticleGrid::Clear()
{
	for (int i = 0; i < n_cells_x; i++)
	{
		for (int j = 0; j < n_cells_y; j++)
		{
			grid[i][j].Clear();
		}
	}
}


NeighbourList::ParticleGrid::~ParticleGrid()
{
	for (int i = 0; i < n_cells_x; i++)
	{
		delete[] grid[i];
	}
	delete[] grid;
}

bool NeighbourList::ParticleGrid::isIndexValid2D(int index_cell_x, int index_cell_y)
{
	return isIndexValid1D(index_cell_x) && index_cell_y >= 0 && index_cell_y < n_cells_y;
}

bool NeighbourList::ParticleGrid::isIndexValid1D(int index_cell_x)
{
	return index_cell_x >= 0 && index_cell_x < n_cells_x;
}

void NeighbourList::ParticleGrid::IndexParticlesOnBoundary(int* index_cell, int size)
{
	if (*index_cell == size)
	{
		*index_cell = size - 1;
	}
}

std::vector<NeighbourList::Particles>  NeighbourList::ParticleGrid::GetNeighbourGridParticles(int index_cell_x, int index_cell_y)
{

	std::vector <Particles> res;

	std::vector<int> result;

	IndexParticlesOnBoundary(&index_cell_x, n_cells_x);
	IndexParticlesOnBoundary(&index_cell_y, n_cells_y);

	auto AddParticlesIfIndexValid
	{
		[&](int x, int y)
		{
			if (isIndexValid2D(x, y))
			{
				res.push_back(grid[x][y]);
			}
		}
	};

	AddParticlesIfIndexValid(index_cell_x - 1, index_cell_y);
	AddParticlesIfIndexValid(index_cell_x, index_cell_y);
	AddParticlesIfIndexValid(index_cell_x + 1, index_cell_y);

	if (parameters::dim == 2)
	{
		AddParticlesIfIndexValid(index_cell_x - 1, index_cell_y - 1);
		AddParticlesIfIndexValid(index_cell_x - 1, index_cell_y + 1);
		AddParticlesIfIndexValid(index_cell_x, index_cell_y - 1);
		AddParticlesIfIndexValid(index_cell_x, index_cell_y + 1);
		AddParticlesIfIndexValid(index_cell_x + 1, index_cell_y - 1);
		AddParticlesIfIndexValid(index_cell_x + 1, index_cell_y + 1);
	}


	return res;
}

double NeighbourList::GetScaleK()
{
	if (parameters::skf == 1)  return 2.0e0;
	else if (parameters::skf == 2 || parameters::skf == 3) return 3.0e0;
	else return 2.E0;
}
double NeighbourList::BlockSize(double hsml)
{
	return GetScaleK() * hsml;
}
int NeighbourList::SizeX(double hsml)
{
	return ceil((parameters::xMaxGeom - parameters::xMinGeomy) / BlockSize(hsml));
}
int NeighbourList::SizeY(double hsml)
{
	return ceil((parameters::yMaxGeom - parameters::yMinGeomy) / BlockSize(hsml));
}
int NeighbourList::Block_x(double *x, double hsml)
{
	double t_index = (x[0] - parameters::xMinGeomy) / BlockSize(hsml);
	return (int)t_index;
}
int NeighbourList::Block_y(double * x, double hsml)
{
	double t_index = (x[1] - parameters::yMinGeomy) / BlockSize(hsml);
	return (int)t_index;
}