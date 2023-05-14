#ifndef HYPERPARAMS_H
#define HYPERPARAMS_H

#define MAX_PARTICLES 1048576

#ifndef MAX_GRID
#define MAX_GRID 1024
#endif
#define MAX_NUM_GRID (MAX_GRID*MAX_GRID)

// #define CUDART_PI_F 3.141592654f

#include <cmath>
#include "MathUtility.hpp"
#include "constants.hpp"

namespace dimensionalConstants
{
  // RodShapedBacterium cell specific constants with dimensions
  constexpr double rodGrwthRtePreFac{ constants::rodGrwthRtePreFac };  // microns / hour
  constexpr double rodModE{ constants::rodModE };                      // Pa - this is currently actually (E/(1-sigma^2))
  constexpr double zeta{ constants::zeta };                      // drag coefficient Pa h
  constexpr double rodSpheroDiam{ constants::rodSpheroDiam };          // microns
  constexpr double avg_div_L{ constants::avg_div_L };            // microns
  constexpr double dt { constants::dt };   // hours
  constexpr double kappa { 1e6 }; // Pa
}
namespace dC = dimensionalConstants;

struct HyperParams
{
  // Particle hypers
  const double radius { 0.5*dC::rodSpheroDiam };
  const double divLength { dC::avg_div_L };
  const double initLength { 0.5*divLength-radius };
  const double growthFactor {
    dC::rodGrwthRtePreFac*dC::zeta / ( dC::rodModE*dC::rodSpheroDiam)
  };

  // Grid hypers
  Uint3 gridSize { 1,1,1 };
  uint numGridCells { 1 };
  Vec3 origin { 0 };
  double max_coll_dist{
    divLength+2*radius
  };
  Vec3 cellSize {
    max_coll_dist,
    max_coll_dist,
    max_coll_dist
  };

  void updateParams(Uint3 _gridSize)
  {
    gridSize = _gridSize;
    numGridCells = gridSize.x*gridSize.y*gridSize.z;
    cellSize = Vec3(
      divLength+2*radius,
      divLength+2*radius,
      divLength+2*radius
    );
    origin = Vec3(
      -0.5*MAX_GRID*cellSize.x,
      -0.5*MAX_GRID*cellSize.y,
      0
    );
  }
};

#endif
