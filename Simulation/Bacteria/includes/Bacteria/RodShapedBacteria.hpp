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

#ifndef ROD_SHAPED_BACTERIA_HPP
#define ROD_SHAPED_BACTERIA_HPP

/*
  RodShapedBacterium particle
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
#include "Geometry.hpp"
// #include "particle.hpp"
#include "IBacterium.hpp"
// #include "Springs.hpp"

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

class RodShapedBacterium : public IBacterium
{

public:
  /*
    Generalised susceptible modelled as a spherocylinder with Hertzian forces
    All quantities are in non-dimensional units
  */

  /*---------------------- Set cell dynamic properties -----------------------*/
  // Cell posn and dynamics
  Vec3 mPos;                          //!< centre of mass vector
  Vec3 mLogPos;                       //!< centre of mass vec at time of binning
  Vec3 mVel{0,0,0};                   //!< centre of mass velocity
  Vec3 mForce{0,0,0};                 //!< net force experienced by this particle
  Vec3 mAngles{0,0.5*constants::pi,0};//!< phi (xy), theta (z), psi (body axis)
  Vec3 mAngVel{0,0,0};                //!< angular velocity in the body frame
  Vec3 mTorque{0,0,0};                //!< net torque experienced by this particle (fixed frame)
  double mLength;                     //!< Length of the cylindircal part

  /*---------------------- Set cell hyperparmeters ---------------------------*/

  static double mAvgGrwthRate;      //!< average growth rate for all cells
  double mGrwthRtePreFac;           //!< cellular growth rate prefactor
  static double mAvgDivLen;         //!< average division length
  static double mRodModE;           //!< proportional to Young's modulus
  static double mRadius;            //!< particle radius
  static uint counter;              //!< counts the total number of bacteria
  uint mId;                         //!< unique id

#if defined(ADHESION)
  static double mKappaDep;          //!< adhesive interaction strength
  static double mRd;                //!< depletion interaction radius
  static double mRi;                //!< repulsion interaction length
  static double mRepStrength;       //!< repulsion interaction strength
#elif defined(CHAINING)
  static double mKappa;             //!< Spring tension
  static double mBendRig;           //!< Bending rigidity
  static double mLinkingProb;       //!< Probability daughters link
  static double mForceThresh;       //!< Threshold force before breaking

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
#elif defined( AG43 )
  static double mKappa;              //!< Spring tension
  static double mForceThresh;        //!< Threshold force before breaking
  // std::vector<IBacterium*> mStuckTo; //!< Springs connect these cells
  SpringHash mSprings; //!< Springs
  virtual SpringHash& getSprings() override
  {
    return mSprings;
  }
#endif


#ifdef PHAGE
  double mTimeSinceInfection{ 0.0 };
  bool mInfectedFlag { false };
  double mLysisPeriod { 1e20 };
  uint mBurstSize { 0 };
  uint mMOI { 0 };

  virtual bool signalLysis() const override
  { return mTimeSinceInfection>=mLysisPeriod; }
  virtual void updateTimeSinceInfection(double dt) override
  { mTimeSinceInfection+=dt; }
  virtual void setInfected() override { mInfectedFlag=true; }
  virtual bool isInfected() override { return mInfectedFlag; }
  virtual void setLysisTime(double _lysis_period) override
  { mLysisPeriod=_lysis_period; }
  virtual double getLysisTime() const override { return mLysisPeriod; }
  virtual void setBurstSize(uint _burst_size) override
  { mBurstSize=_burst_size; }
  virtual uint getBurstSize() const override { return mBurstSize; }
  virtual void updateMOI() override { ++mMOI; }
  virtual uint getMOI() const override { return mMOI; }
#endif

  uint getID() const override
  {
    return mId;
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
    // const double effective_length_3{ mRadius*mRadius*( mRadius + 0.75*mLength ) };
    // const double effective_length{ pow(effective_length_3,(1.0/3.0)) };
    // mVel=mForce/effective_length;
    mVel=mForce/mLength;
  }

  // Return the lab frame angvel
  virtual Vec3& getAngVel() override
  {
    return mAngVel;
  }

  // Set the lab frame angvel from the torque on this rod
  virtual void setAngVel() override
  {
#ifndef ANISOTROPIC
    // const double effective_length_3{ mRadius*mRadius*( mRadius + 0.75*mLength ) };
    // mAngVel=(4.0/3.0)*mTorque/effective_length_3;
     mAngVel = mTorque * 12 / (mLength*mLength*mLength);
#else
    std::cout << "exit from set angvel rod shaped" << '\n';
    exit(1);
#endif
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
  RodShapedBacterium (
    double rcm_x=0,
    double rcm_y=0,
    double rcm_z=0,
    double theta=0,
    double alpha=constants::pi*0.5,
    double grwthPreFac=constants::nondim_rodGrwthRtePreFac,
    double init_length=constants::nondim_init_length
  );
  RodShapedBacterium (
    const Vec3 &rcm,
    double theta=0,
    double alpha=constants::pi*0.5,
    double grwthPreFac=constants::nondim_rodGrwthRtePreFac,
    double init_length=constants::nondim_init_length
  );

  /*--------------------- --------------------------*/
  virtual std::string getMyType() const override
  {
#ifdef CHAINING
    return "ChainingRodShaped";
#elif defined(AG43)
  return "AG43RodShaped";
#else
    return "RodShaped";
#endif
  }

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
    mAngles.x=theta;
    mAngles.y=alpha;
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
    return mRodModE;
  }

  /*--------------------- Specific cell abilities ----------------------------*/

  // double consumeNutrients(double c = constants::conc_max,
  //                         double c_half = constants::c_sat);
  /**<
  \brief Find the rate bacteria consume the local nutrients
  @param[in] c: local concentration of nutrients
  @param[in] c_half: half-saturation of the Monod function
  @returns rate of nutrient consumption
  */

  virtual void grow(double dt) override;
  /**<
  \brief Grow the bacteria based on a constant nutrient density. This is grows bacteria
  as in [2];
  @param[in] dt: simulation timestep
  @param[out] mLength of the bacteria is increased according to growth rate
  @returns Void.
  */

  // void grow(double c = constants::conc_max, double dt = 0.5e-4);
  /*
  Grow the bacteria based on the local nutrient availability. This has the effect
  of increasing mLength as in [1]. This expression also accounts for the
  local nutrient uptake of different sizes of bacteria.
  Parameters:
    c: local concentration of nutrients
    dt: simulation timestep
    Effect:
      Updates mLength of the bacteria
  */

  virtual bool signalDivide() override
  {
    return mLength>mAvgDivLen;
  }
  /**<
  \brief Check if cell will divide.
  @returns True (if this is the case)
  */

  virtual void divide(std::vector<IBacterium*>& cell_list) override;

  // double getDaughterTheta();
  /**<
  \brief Get the orientation of the daughter cell in the plane
  @returns theta: mother theta plus a small kick due to fluctations
  */

  // double getDaughterAlpha();
  /**<
  \brief Get the orientation of the daughter cell along z
  @returns alpha: mother alpha plus a small kick due to fluctations
  */

  // double getDaughterGrowthRate();
  /**<
  \brief Get the growth rate of the daughter cell
  @returns grwth_rate: uniformly distributed growth rate in the range 0.5-1.5 of average
  */

  virtual void move(double dt) override;
  /**<
  \brief Move the RodShapedBacterium according to Newtonian dynamics.
  @param[in] dt: timestep
  @param[out] mPos: translated according to forces on this particle
  @param[out] theta,alpha: rotated according to torques on this particle
  */

  // void setRCMVelFromForce();
  // void setRCMAngVelFromTorque();

  void reset() override;      // Zero any dynamical quantity

  double getCellArea();
  /*
  Find the cell area based on spherocylinder shape
  Returns :
    cell area
  */

  virtual Vec3 getOrientation() const override;
  /**<
    \brief Find vector of orientation.
    @returns (sin(alpha) cos(theta), sin(alpha) sin(theta), cos(alpha) )
  */

  virtual void getMyEndVecs(Vec3& p,Vec3& q) const override;
  /**<
    \brief Find the two end point vectors of the spherocylinder
    vec(P1,P2): the start and end point vectors respectively
    @returns void
  */

  /**
    \brief Find radius of the enclosing sphere this particle
    @return effective_radius
  */
  virtual double getEffectiveR() const override
  {
    return mRadius+0.5*mLength;
  }

  virtual double getRadius() const override
  {
    return mRadius;
  }

  virtual double getLength() const override
  {
    return mLength;
  }

  virtual ~RodShapedBacterium () {};

  // ----------------------------------- IO ------------------------------------
  virtual void printToFile(std::ostream &out) override
  {
#ifdef CHAINING
    out << getMyType()       << "\t";
    out << mId               << "\t";
    out << mLength           << "\t";
    out << mRadius           << "\t";
    out << mPos              << "\t";
    out << getOrientation()  << "\t";
    if (getLowerLink()) out << getLowerLink()->getID() << "\t";
    else out << "None" << "\t";
    if (getUpperLink()) out << getUpperLink()->getID() << "\n";
    else out << "None" << "\n";
#elif defined(AG43)
    out << getMyType()       << "\t";
    out << mId               << "\t";
    out << mLength           << "\t";
    out << mRadius           << "\t";
    out << mPos              << "\t";
    out << getOrientation()  << "\t";
    out << "[";
    for ( auto &spring : mSprings )
    {
      const Springs& ss { spring.second };
      assert( ss.mCellB->getID()==spring.first );
      out << "("
          << ss.mCellB->getID() << "," << ss.mS << "," << ss.mT << ","
          << ss.mOriLink.x << "," << ss.mOriLink.y << "," << ss.mOriLink.z
          << "),";
    }
    out << "]\n";
#else
    out << getMyType()       << "\t";
    out << mId               << "\t";
    out << mLength           << "\t";
    out << mRadius           << "\t";
    out << mPos              << "\t";
    out << getOrientation()  << "\n";
#endif
  }
};

// double getDaughterAngle(double angle);
/**<
  slightly perturb input angle to give quasi-spherical shape
  Parameters:
    angle: the mother's angle
  Returns:
    daughters perturbed angle
*/


// void createDaughters(
//   uint ii,
//   std::vector<RodShapedBacterium> &cell_list
// );

// RodShapedBacterium::cellpair createDaughters(RodShapedBacterium &mother_cell);
// void createDaughters(RodShapedBacterium &mother_cell,
//                      std::vector<RodShapedBacterium> &daughter_cells);
/**<
 Create two daughter cells with slightly perturbed orientation and growth rate
 compared to the mother. The daughters' centres are placed at rcm +/- 0.25*(L+d)
 with length (L-d)/2
 Effect:
   Emplaces two daughter cells at the end of the passed daughter_cells vector
*/

void getEndVecs(const RodShapedBacterium &cell, Vec3 &p, Vec3 &q);
/**<
  Find the two end point vectors of a RodShapedBacterium
  Returns:
    vec(P1,P2): the top and bottom end point vectors respectively
*/

std::ostream& operator<< (std::ostream &out, const RodShapedBacterium &cell);
/**<
  Create representation of cell for easy I/O
  This will print in the format
  rcm (rx,ry,rz) orientation (ox,oy,oz) length L
*/

std::ostream& printCellDataToFile(std::ostream &out,
                                  const RodShapedBacterium &cell);

void printInteraction(std::ostream &out, const RodShapedBacterium &cell,
                      const Vec3 &virt_cen,
                      const Vec3 &force_on_me,
                      const Vec3 &torque_on_me);

void initialiseRodParameters(double aspect_ratio, double growth_rate);
/**<
 * Set the hyper parameters for the general rod shaped bacterium class
 */

#ifdef ADHESION
void initialiseAdhesionParameters(
  double kappa_depletion,
  double Rd,
  double Ri,
  double rep_strength
);
/**<
  Check the hyperameters for adhesion bacteria are suitable and set the
  associated static variables.
*/
#endif

#ifdef CHAINING
void initialiseChainingParameters(
  double kappa,
  double bend_rig,
  double linking_prob,
  double force_thresh
);
/**<
  Check the hyperameters for chaining bacteria are suitable and set the
  associated static variables.
*/
#endif

#ifdef AG43
void createSpringLinks(
  std::vector<IBacterium*> &cells
);
/**<
  @brief Set up spring links between the bacteria if they come into contact.
*/

void removeSpringLinks(
  std::vector<IBacterium*> &cells
);
/**<
  @brief Remove spring links between the bacteria if they get too far away.
*/

void initialiseAg43Parameters(
  double div_len,
  double growth_rate,
  double kappa,
  double force_thresh
);
/**<
  @brief Initialise ag43 hyperparmeters
*/

void checkSpringLinks(
  std::vector<IBacterium*> &cells
);
/**<
  @brief Check all springs are connected between the bacteria in the right way
*/
#endif

#endif // End fileguard
