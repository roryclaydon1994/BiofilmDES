// Standard libraries
#include <iostream>
#include <cmath>
#include <vector>
#include <array>
#include <random>
#include <memory>

// User defined libraries
#include "constants.hpp" // definition of constants namespace
#include "MathUtility.hpp"
#include "particle.hpp"

/*-------------------------- Representations ---------------------------------*/
std::ostream& operator<< (std::ostream &out, const Particle &particle)
{
  out << particle.mId     << "\t";
  out << particle.mPos;
  return out;
}

// Define static variables
std::mt19937 Particle::mGenerator { constants::SEED };
// long Particle::counter { 0 }; //!< unique particle counter
std::array<std::string,4> Particle::col_headers
{
  "particle_id\t",
  "com_vec_x\t",
  "com_vec_y\t",
  "com_vec_z"
};
