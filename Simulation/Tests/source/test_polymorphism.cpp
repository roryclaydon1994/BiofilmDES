/*******************************************************************************

Script to test polymorphism

*******************************************************************************/
// Standard libraries
#include <iostream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <array>
#include <random>
#include <memory>
#include <cassert>
#include <fstream>

// User defined libraries
#include "constants.hpp"                    // definition of constants namespace
#include "MathUtility.hpp"                   // Vec3 defined here
#include "particle.hpp"                     // Trial a new particle class
#include "VerletGrid.hpp"
#include "VerletGrid.cpp"

template class GridCell<Particle>;

class FakeParticleA : public Particle
{
public:
  FakeParticleA (double rcm_x=0, double rcm_y=0, double rcm_z=0) :
  Particle {rcm_x,rcm_y,rcm_z}
  {}

  ~FakeParticleA ()
  {}
};


class FakeParticleB : public FakeParticleA
{
public:
  FakeParticleB (double rcm_x=0, double rcm_y=0, double rcm_z=0) :
  FakeParticleA {rcm_x,rcm_y,rcm_z}
  {}

  ~FakeParticleB ()
  {}
};

void testBasics ()
{
  Particle MrP;
  FakeParticleA fpA;
  FakeParticleB fpB;

  std::vector< Particle* > poly_test_vec {&MrP,&fpA,&fpB};

  for ( auto &pp : poly_test_vec )
  {
    std::cout << *pp << '\n';
  }
}

void testGrid() {

  // Create some random particles
  std::vector< Particle > particles;
  std::vector< FakeParticleA > fakeAs;
  std::vector< FakeParticleB > fakeBs;
  std::normal_distribution<double> distribution(0.0, 1.0);

  for ( int ii=0; ii<10; ++ii )
  {
    particles.emplace_back(
      distribution(Particle::mGenerator),
      distribution(Particle::mGenerator),
      std::abs(distribution(Particle::mGenerator))
    );

    fakeAs.emplace_back(
      distribution(FakeParticleA::mGenerator),
      distribution(FakeParticleA::mGenerator),
      std::abs(distribution(FakeParticleA::mGenerator))
    );

    fakeBs.emplace_back(
      distribution(FakeParticleB::mGenerator),
      distribution(FakeParticleB::mGenerator),
      std::abs(distribution(FakeParticleB::mGenerator))
    );
  }

  // Create a vector to control all particles
  std::vector< Particle* > poly_test_vec;
  poly_test_vec.reserve(30);

  for ( auto &pp : particles )
  {
    poly_test_vec.push_back(&pp);
  }
  for ( auto &pp : fakeAs )
  {
    poly_test_vec.push_back(&pp);
  }
  for ( auto &pp : fakeBs )
  {
    poly_test_vec.push_back(&pp);
  }


  // Boxes must have at least a side length of the maximum distance between cells
  GridCell<Particle>::box_width = 1;
  std::vector< std::unique_ptr<GridCell<Particle>> > grid_cells;


  std::cout << "1" << '\n';
  allocateGrid(grid_cells,poly_test_vec);
  std::cout << "2" << '\n';
  populateNeighbourCells(grid_cells);        // Populate grid neighbour list
  std::cout << "3" << '\n';
  updateGridCellLists(grid_cells,poly_test_vec); // Bin cells to grid
  std::cout << "4" << '\n';
  createNeighbourLists(grid_cells,poly_test_vec);// Create cell neighbour lists
  std::cout << "5" << '\n';

  for ( auto &pp : fakeAs )
  {
    std::cout << "---------------------------------------------------" << '\n';
    std::cout << pp.mId << '\n';
    for ( auto &nn : pp.mMyNeighbourList )
    {
      std::cout << *nn << '\n';
    }
  }
}

int main(int argc, char const *argv[]) {
  // testBasics();
  testGrid();
  return 0;
}
