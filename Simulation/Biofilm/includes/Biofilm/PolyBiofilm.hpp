#ifndef POLYBIOFILM_HPP
#define POLYBIOFILM_HPP

// Standard libraries
#include <vector>

// Custom
#include "constants.hpp"
#include "IBacterium.hpp"
#include "VerletGrid.hpp"
#include "Phage.hpp"
#include "RodShapedBacteria.hpp"

class PolyBiofilm
{
public:
  /*----------------------------- Sim Parameters -----------------------------*/
  double mDt { constants::nondim_dt };
  std::vector<IBacterium*> mCells;
#ifdef PHAGE
  std::vector<Phage*> mPhage;
  uint mNumberOfInfected { 0 };
#endif
  uint  mOutFreq;
  double mTargetSize{ constants::colony_size_target }; // If length > this, stop

  ulong mTimeSteps;
  GridCells<IBacterium> mGrid;

  /*------------------------------- Utilities --------------------------------*/
  PolyBiofilm (
#ifndef PHAGE
    std::vector<IBacterium*>& _cells,
#else
    std::vector<IBacterium*>& _cells,
    // std::vector<Phage*>& _phage
#endif
    double _box_width=constants::box_width
  ) :
    mDt { constants::nondim_dt },
    mCells{ _cells },
#ifdef PHAGE
    // mPhage{ _phage },
#endif
    mGrid{ _cells, _box_width }
  {

    std::cout << "The target size is " << mTargetSize << '\n';

    std::cout << "box size\t: " << mGrid.mCellSize << " vs "
              << constants::box_width << '\n';

    // Number of timesteps between each output
    // Previously 0.02
    mOutFreq   = static_cast<uint>( ceil( ( 0.1/constants::baseTimeScale )/mDt) );

    // Number of simulation hours to run for
    mTimeSteps = static_cast<ulong>( ceil( ( 24/constants::baseTimeScale )/mDt) );

    if ( mOutFreq > mTimeSteps )
    {
      std::cout << "Error - output step size less than time_steps!" << '\n';
      std::cout << "Output frequency: " << mOutFreq
                << " timesteps: " << mTimeSteps
                << '\n';
      std::cout << "Exiting..." << '\n';
      exit (EXIT_FAILURE);
    }

    std::cout << "outputting every " << mOutFreq*mDt
              << " or every " << mOutFreq*mDt*constants::baseTimeScale
              << " hours "  << '\n';
  }

  void createLogFile();

  void updateStatus(long output_counter,ulong num_outputs)
  {
    std::cout << "\rsaving at time: "
              << std::left << std::setw(4) << std::setprecision(3)
              << output_counter*mOutFreq*mDt*constants::baseTimeScale
              << " hrs output: "
              << std::setw(5) << std::right << output_counter
              << " / "
              << std::left << std::setw(5) << num_outputs
              << " ---- "
              << "num bacteria: " << std::setw(6) << std::left << mCells.size()
#ifdef PHAGE
              << " "
              << "num infected: " << std::setw(6) << std::left << mNumberOfInfected
              << " "
              << "num phage: "    << std::setw(6) << std::left << mPhage.size()
#endif
              << std::flush;
  }

  void updateOneTimeStep(bool& update_neighbours,uint &verlet_counter);

  void runSim();

  /**
    Get the total "length" of the colony
  */
  double getLength()
  {
    double length{ 0.0 };
    for ( auto cell : mCells )
      length+=cell->getLength()+2*cell->getRadius();
    return length;
  }

  ~PolyBiofilm()
  {
    for ( auto &cell : mCells )
    {
      delete cell;
    }
#ifdef PHAGE
  for ( auto &phage : mPhage ) delete phage;
#endif
  }
};

#endif
