// Standard libraries
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <cassert>
#include <random>
#include <memory>

// User defined
#include "SphericalBacteria.hpp"
#include "MathUtility.hpp"
#include "IO.hpp"
#include "constants.hpp"        // definition of constants namespace

double SphericalBacterium::mAvgDivRad{
  constants::nondim_avg_div_radius
};                                                  // average division radius
double SphericalBacterium::mModE
{
  constants::nondim_sphericalModE
};                                                  // Young's modulus

uint SphericalBacterium::counter { 0 };             // unique counter
