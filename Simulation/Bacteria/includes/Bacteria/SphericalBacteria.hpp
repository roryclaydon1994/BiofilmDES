/*******************************************************************************
    BiofilmDES  - a program that simulates a growing colony of microbial cells

    Contributing author:
    Rory Claydon, University of Edinburgh, rory.claydon@ed.ac.uk

    Copyright (2020) The University of Edinburgh.

    The software is based on algorithms described in:

    Mechanically driven growth of quasi-two dimensional microbial colonies,
    F.D.C. Farrell, O. Hallatschek, D. Marenduzzo, B. Waclaw,
    Phys. Rev. Lett. 111, 168101 (2013).

    Three-dimensional distinct element simulation of spherocylinder crystallization.
    Pournin, L., Weber, M., Tsukahara, M. et al.
    Granul. Matter 7, 119–126 (2005).

    A fast algorithm to evaluate the shortest distance between rods,
    C. Vega, S. Lago,
    Comput. Chem., 18(1), 55-59 (1994)

    I would like to thank Bartlomiej Waclaw from Edinburgh University for some
    very useful discussions on algorithm stability, timestep choice and some
    potential optimisations to try out in future.

    This file is part of BiofilmDES.

    BiofilmDES is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    BiofilmDES is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
    See the GNU General Public License for more details.

    A copy of the GNU General Public License can be found in the file
    License.txt or at <http://www.gnu.org/licenses/>.

    Compilation and run from current directory:
      make && ./biofilm.out 0 1.1 0.95

    Further details in the documentation

*******************************************************************************/

#ifndef SPHERICAL_BACTERIA_HPP
#define SPHERICAL_BACTERIA_HPP

/*
  SphericalBacterium particle
*/

// Standard libraries
#include <iostream>
#include <cmath>
#include <vector>
#include <array>
#include <random>
#include <memory>
#include <unordered_map>

// User defined libraries
#include "constants.hpp" // definition of constants namespace
#include "MathUtility.hpp"
#include "IBacterium.hpp"

/*----------------------------------------------------------------------------*/
/*
  Bacteria base class defined here for use in agent based rodModElling of biofilms
  The parameter values and simulaton method for the biofilm component of the
  code is strongly based on the following article:
  [1] O. Hallatschek et al., Phys. Rev. Lett., 111, 168101 (2013)
      DOI:https://doi-org.ezproxy.is.ed.ac.uk/10.1103/PhysRevLett.111.168101
  [2] Luca Giomi
  Author: Rory Claydon
*/
/*----------------------------------------------------------------------------*/

class SphericalBacterium : public IBacterium
{

public:

  /*---------------------- Set cell dynamic properties -----------------------*/
  Vec3 mPos;                          //!< centre of mass vector
  Vec3 mLogPos;                       //!< centre of mass vec at time of binning
  Vec3 mVel{0,0,0};                   //!< centre of mass velocity
  Vec3 mForce{0,0,0};                 //!< net force experienced by this particle
  Vec3 mAngles{0,0,0};                //!< phi (xy), theta (z), psi (body axis)
  Vec3 mAngVel{0,0,0};                //!< angular velocity (body frame)
  Vec3 mTorque{0,0,0};                //!< net torque experienced by this particle

  /*---------------------- Set cell hyperparmeters ---------------------------*/

  double mGrwthRtePreFac;                //!< cellular growth rate prefactor
  static double mAvgDivRad;              //!< average division radius
  static double mModE;                   //!< The Young's modulus
  static uint counter;                   //!< counts the total number of bacteria
  double mRadius{
    constants::nondim_spherical_init_radius
  };                                     //!< Cell radius
  uint mId;                              //!< unique id

  uint getID() const override
  {
    return mId;
  }

  /*--------------------- --------------------------*/

#ifdef CHAINING
  // Pointers to who this cell is connected to
  IBacterium* mUpperEndLinkedTo { nullptr };
  IBacterium* mLowerEndLinkedTo { nullptr };
  virtual void setUpperLink(IBacterium* cell) override
  { mUpperEndLinkedTo=cell; }
  virtual void setLowerLink(IBacterium* cell) override
  { mLowerEndLinkedTo=cell; }
  virtual IBacterium* getUpperLink() const override
  { return mUpperEndLinkedTo; }
  virtual IBacterium* getLowerLink() const override
  { return mLowerEndLinkedTo; }
#endif

  virtual Vec3 getPos() const override
  {
    return mPos;
  }

  virtual void setPos(double x, double y, double z=0.0) override
  {
    mPos.x=x; mPos.y=y; mPos.z=z;
  }

  virtual void setAngles(double theta,double alpha=0.5*constants::pi) override
  {
    std::cout << "Error! Do not set angles for spherical bacteria" << '\n';
    exit(67);
  }

  virtual Vec3 getLoggedPos() const override
  {
    return mLogPos;
  }

  virtual void setLoggedPos() override
  {
    mLogPos=mPos;
  }

  virtual double getModE() const override
  {
    return mModE;
  }

  virtual Vec3& getForce() override
  {
    return mForce;
  }
  virtual void addForce(Vec3 force) override
  {
    mForce+=force;
  }
  virtual Vec3& getTorque() override
  {
    return mTorque;
  }
  virtual void addTorque(Vec3 torque) override
  {
    mTorque+=torque;
  }
  virtual Vec3& getVel() override
  {
    return mVel;
  }
  virtual void setVel() override
  {
    mVel=mForce/mRadius;
  }
  virtual Vec3& getAngVel() override
  {
    return mAngVel;
  }
  virtual void setAngVel() override
  {
    mAngVel=(4.0/3.0)*mTorque/pow(mRadius,3);
  }

  /* -------------------------- Cell Linked List -----------------------------*/
  uint mMyGridCell { constants::INF };
  IBacterium* mHeadLink { nullptr };
  std::vector<IBacterium*> mMyNeighbourList; //!< neighbour list
  virtual std::vector<IBacterium*>& getNeighbourList()
  {
    return mMyNeighbourList;
  }
  virtual IBacterium* getHeadLink() override
  {
    return mHeadLink;
  }
  virtual void setHeadLink(IBacterium* new_head) override
  {
    mHeadLink = new_head;
  }
  virtual uint& getGridCell()
  {
    return mMyGridCell;
  }


  /*--------------------------- Constructors ---------------------------------*/
  SphericalBacterium (
    double rcm_x=0,
    double rcm_y=0,
    double rcm_z=0,
    double grwthPreFac=constants::nondim_sphericalGrwthRtePreFac,
    double radius=constants::nondim_spherical_init_radius
  ) : mPos{rcm_x,rcm_y,rcm_z}, mGrwthRtePreFac{grwthPreFac},mRadius{radius}
  {
    mId = counter++;           // increment unique counter
  }

  SphericalBacterium (
    const Vec3 &rcm,
    double grwthPreFac=constants::nondim_sphericalGrwthRtePreFac,
    double radius=constants::nondim_spherical_init_radius
  ) : SphericalBacterium {rcm.x,rcm.y,rcm.z,grwthPreFac,radius} {}

  virtual std::string getMyType() const override
  {
    return "Spherical";
  }

  void grow(double dt) override
  {
    mRadius+=dt*mGrwthRtePreFac;
  }
  /**<
  \brief Grow the bacteria based on a constant nutrient density. This is grows bacteria
  as in [2];
  @param[in] dt: simulation timestep
  @param[out] mLength of the bacteria is increased according to growth rate
  @returns Void.
  */

  virtual bool signalDivide() override
  {
    return mRadius>mAvgDivRad;
  }

  virtual void divide(std::vector<IBacterium*>& cell_list) override;

  virtual void move(double dt) override
  {
#ifndef MOVE_3D
    mVel.z=0;
#endif
    mPos+=dt*mVel;
  };
  /**<
  \brief Move the SphericalBacterium according to Newtonian dynamics.
  @param[in] dt: timestep
  @param[out] mPos: translated according to forces on this particle
  @param[out] theta,alpha: rotated according to torques on this particle
  */

  void reset() override
  {
    mVel.zero();
    mAngVel.zero();
    mForce.zero();
    mTorque.zero();
  }

  Vec3 getOrientation() const override
  {
    std::cout << "Error: SphericalBacterium has no orientation" << '\n';
    exit (EXIT_FAILURE);
  }
  /**<
    \brief Find vector of orientation.
    @returns (sin(alpha) cos(theta), sin(alpha) sin(theta), cos(alpha) )
  */

  void getMyEndVecs(Vec3 &p, Vec3 &q) const override
  {
    p=q=mPos;
  }
  /**<
    \brief Find the two end point vectors of the spherocylinder
    @returns vec(P1,P2): the top and bottom end point vectors respectively
  */

  /**
    \brief Find radius of the enclosing sphere this particle
    @return effective_radius
  */
  double getEffectiveR() const override
  {
    return mRadius;
  }

  virtual double getRadius() const override
  {
    return mRadius;
  }
  virtual double getLength() const override
  {
    return 0.0;
  }

  virtual ~SphericalBacterium () {};

  // ----------------------------------- IO ------------------------------------
  virtual void printToFile(std::ostream &out) override
  {
    out << getMyType()       << "\t";
    out << mId               << "\t";
    out << 0                 << "\t";
    out << mRadius           << "\t";
    out << mPos              << "\t";
    out << Vec3{0.0,0.0,0.0} << "\n";
  }
};

#endif // End fileguard
