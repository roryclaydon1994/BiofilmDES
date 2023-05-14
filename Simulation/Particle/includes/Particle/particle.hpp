
#ifndef PARTICLE
#define PARTICLE

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

class Particle
{
public:
  /*-------------------- Set particle dynamic properties ---------------------*/
  // Cell posn and dynamics
  Vec3 mPos;                        //!< centre of mass vector
  Vec3 mVel{0,0,0};              //!< centre of mass velocity
  Vec3 mForce{0,0,0};      //!< net force experienced by this particle

  /*----------------------- Set particle hyperparmeters ----------------------*/

  // Create mersenne twister generator to control random numbers associated with
  // instances of this class
  static std::mt19937 mGenerator;     //!< All instances share same generator

  // static long counter;       //!< counts the total number of susceptibles
  long mId { 0 };               //!< unique id of particle

  static std::array<std::string,4> col_headers;

  /*-------------------- Cell linked list properties -------------------------*/

  //! Pointer to previous particle in list - never do resource management with this
  Particle* mHeadLink { nullptr };
  std::vector<Particle*> mMyNeighbourList; //!< neighbour list

  //! The index of the grid this particle lives in
  long mMyGridCell { -1 };

  /**
    \brief Find radius of the enclosing sphere this particle
    @return effective_radius
  */
  virtual double getEffectiveR()
  {
    return 0.5*constants::nondim_rodSpheroDiam;
  }

  /*-------------------- Constructors -------------------------*/
  Particle (double rcm_x=0, double rcm_y=0, double rcm_z=0) :
    mPos {rcm_x,rcm_y,rcm_z}
  {
    // mId = counter++;           // increment unique counter
  }

  Particle (const Vec3 &rcm) :
    Particle(rcm.x,rcm.y,rcm.z)
  {}

  ~Particle ()
  {}

  /*----------------------- Dynamic Functions --------------------------------*/
  virtual void move(double dt)
  {
#ifdef CONFINE_2D
    mVel.z = 0.0;
#endif // End 2D
    mPos += dt*mVel;
  }

  virtual void reset()
  {
    mVel.x         =  mVel.y         = mVel.z         = 0;
    mForce.x =  mForce.y = mForce.z = 0;
  }

};

std::ostream& operator<< (std::ostream &out, const Particle &particle);

#endif
