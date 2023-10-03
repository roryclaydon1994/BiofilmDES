// Standard libraries
#include <iostream>
#include <cmath>
#include <vector>
#include <array>
#include <cassert>
#include <random>
#include <chrono>
#include <memory>
#include <numeric>

// Custom classes
#include "constants.hpp"         // definition of constants namespace
#include "MathUtility.hpp"
#include "IBacterium.hpp"
#include "RodShapedBacteria.hpp"
#include "SphericalBacteria.hpp"
#include "VerletGrid.hpp"

std::mt19937 generator {
  static_cast<unsigned long>(
    std::chrono::high_resolution_clock::now().time_since_epoch().count()
  )
};

class Pseudomonas : public RodShapedBacterium
{
public:
  Pseudomonas (Vec3 _rcm) :
  RodShapedBacterium {_rcm}
  {}

  virtual std::string getMyType() const override
  {
    return "Pseudomonas";
  }

  ~Pseudomonas ()
  {}

};

template < class S >
void setUp(std::vector<S*> &sus)
{
  std::normal_distribution<double> distribution(-10.0, 10.0);
  for ( int ii=0; ii<100; ++ii )
  {
    sus.push_back(
      new RodShapedBacterium(
        Vec3(
          distribution(generator),
          distribution(generator),
          distribution(generator)
        )
      )
    );
    sus.push_back(
      new Pseudomonas(
        Vec3(
          distribution(generator),
          distribution(generator),
          distribution(generator)
        )
      )
    );
    sus.push_back(
      new SphericalBacterium(
        Vec3(
          distribution(generator),
          distribution(generator),
          distribution(generator)
        )
      )
    );
  }
}

void testGrid(std::vector<IBacterium*> &sus)
{
  GridCells grid{ sus, constants::nondim_rodSpheroDiam + constants::nondim_avg_div_L };
  grid.updateGridCellLists(sus);
  grid.createNeighbourLists(sus);

  for ( auto &particle : sus )
  {
    if ( particle->getNeighbourList().size() > 0 )
    {
      std::cout << "P: " << particle->getMyType() << '\n';
      for ( auto &neighbour_particle : particle->getNeighbourList() )
      {
        std::cout << "N:" << neighbour_particle->getID()
                  << " " << neighbour_particle->getMyType() << '\n';
      }
    }
  }
}

int main(int argc, char const *argv[])
{
  std::vector< IBacterium* > sus;
  setUp(sus);
  testGrid(sus);
  return 0;
}
