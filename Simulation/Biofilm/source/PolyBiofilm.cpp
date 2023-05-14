// Standard libraries
#include <iostream>
#include <cmath>
#include <vector>
#include <array>
#include <cassert>
#include <random>
#include <chrono>
#include <string>
#include <memory>
#include <numeric>
#include <filesystem>
// //#include <execution>

// Custom classes
#include "constants.hpp"         // definition of constants namespace
#include "MathUtility.hpp"
#include "IBacterium.hpp"
#include "RodShapedBacteria.hpp"
#include "RandUtil.hpp"
#include "VerletGrid.hpp"
#include "forces.hpp"
#include "IO.hpp"
#include "PolyBiofilm.hpp"

void PolyBiofilm::createLogFile()
{
  std::filesystem::create_directories(sim_out_dir);
  std::string filename=sim_out_dir+'/'+"biofilm.log";
  std::ofstream out_file;
  out_file.open( filename, std::ios::out );
  if (!out_file)
  {
    std::cerr << "Failed to open "
              << filename
              << " for writing!"
              << std::endl;
    exit (10);
  }
  std::vector<std::string> headers{};
  std::vector<double> values{};

  /*=== save basics ===*/
  headers.push_back("OutputFrequency");
  values.push_back(mOutFreq);

  headers.push_back("RodAspectRatio");
  values.push_back(RodShapedBacterium::mAvgDivLen/(2*RodShapedBacterium::mRadius));

  headers.push_back("RodModE");
  values.push_back(RodShapedBacterium::mRodModE);

  // In microns / hr
  headers.push_back("RodGrowthRate");
  values.push_back(RodShapedBacterium::mAvgGrwthRate);

#if defined(ADHESION)
  headers.push_back("KappaDep");
  values.push_back(RodShapedBacterium::mKappaDep);

  headers.push_back("Rd");
  values.push_back(RodShapedBacterium::mRd);

  headers.push_back("Ri");
  values.push_back(RodShapedBacterium::mRi);

  headers.push_back("RepStrength");
  values.push_back(RodShapedBacterium::mRepStrength);
#elif defined(CHAINING)
  headers.push_back("Kappa");
  values.push_back(RodShapedBacterium::mKappa);

  headers.push_back("BendRig");
  values.push_back(RodShapedBacterium::mBendRig);

  headers.push_back("LinkingProb");
  values.push_back(RodShapedBacterium::mLinkingProb);

  headers.push_back("ForceThresh");
  values.push_back(RodShapedBacterium::mForceThresh);
#elif defined(AG43)
  headers.push_back("Kappa");
  values.push_back(RodShapedBacterium::mKappa);

  headers.push_back("ForceThresh");
  values.push_back(RodShapedBacterium::mForceThresh);
#endif

  for ( uint ii=0; ii<headers.size(); ++ii )
  {
    out_file << headers[ii] << '\t';
  }
  out_file << "SEED" << '\n';

  for ( uint ii=0; ii<values.size(); ++ii )
  {
    out_file << values[ii] << '\t';
  }
  out_file << gen_rand.getSeed() << '\n';
  out_file.close();
  std::cout << "Saved log file to " << filename << '\n';
}

void PolyBiofilm::updateOneTimeStep(bool &update_neighbours, uint &verlet_counter)
{
  // ---------------------- Grow and divide cells ------------------------------
  // std::cout << "tries to divide " << mCells.size() << '\n';
  for ( int ii=mCells.size()-1; ii>=0; --ii )
  {
#ifdef PHAGE
    if ( mCells[ii]->isInfected()==false )
#endif
    {
      mCells[ii]->grow(mDt);
      if ( mCells[ii]->signalDivide() )
      {
        // std::cout << "attempts division" << '\n';
        mCells[ii]->divide(mCells);
        // std::cout << "completes division" << '\n';

        // Any division will require a re-binning of the cells
        update_neighbours=true;
        // std::cout << "Binning: division" << '\n';
      }
    }
  }

  // ---------------------- Generate Verlet List ------------------------------
  // cell lists need to be updated
  if ( update_neighbours )
  {
    // std::cout << "updating nieghbour list " << mCells.size() << '\n';
    mGrid.updateVerletList(mCells);
    verlet_counter=0;
    update_neighbours=false;
  }

#ifdef AG43
  // Create springs between the cells if they come into contact
  createSpringLinks(mCells);
#endif

  // std::cout << "interacts cells: " << mCells.size() << '\n';
  // Update the force and torque for each cell to the current time
  polyInteractParticles(mCells);

#ifdef AG43
  // Remove links if they get too far away
  removeSpringLinks(mCells);
#endif

  // std::cout << "enters update" << '\n';
  // std::cout << "move and grow" << '\n';
  for ( auto &cell : mCells )
  {
    // std::cout << cell->getID()     << " "
    //           << cell->getMyType() << " "
    //           << cell->getPos()    << " "
    //           << cell->getRadius()
    //           << '\n';
    cell->move(mDt);
    double moved {
      dot2( cell->getLoggedPos()-cell->getPos() )
    };
    if ( moved>dot2( mGrid.mVerletMoveR ) )
    {
      // std::cout << "Binning: moved too far" << '\n';
      update_neighbours=true;
    }
    cell->reset();
  }

#ifdef PHAGE
  // Add the phage interaction here
  // phage at their current positions attempt to infect the bacteria
  // move remining phage
  // std::cout << "\n 1 phage infections" << '\n';
  attemptInfections(mCells,mPhage,mGrid);
  // std::cout << "\n 2 phage infections" << '\n';

  // std::cout << "\n 1 move phage" << '\n';
  for ( auto phage : mPhage ) { phage->move(mDt); }
  // std::cout << "\n 2 move phage" << '\n';

  // std::cout << "\n 1 lyse cells" << '\n';
  mNumberOfInfected=0;
  for ( int ii=mCells.size()-1; ii>=0; --ii )
  {
    if ( mCells[ii]->isInfected() )
    {
      ++mNumberOfInfected;
      mCells[ii]->updateTimeSinceInfection(mDt);
      if ( mCells[ii]->signalLysis() )
      {
        // std::cout << "lysing cell " << mCells.size() << '\n';
        lyseCell(mCells[ii],mPhage);   // Create phage uniformly over lysed cell
        delete mCells[ii];
        mCells.erase(mCells.begin()+ii);  // Remove the dead cell
        // std::cout << "now " << mCells.size() << '\n';
        update_neighbours=true;
      }
    }
  }
  // std::cout << "\n 2 lyse cells" << '\n';
#endif

#if defined(AG43) && !defined(NDEBUG)
  checkSpringLinks(mCells);
#endif

}

void PolyBiofilm::runSim()
{
  createLogFile();

  std::cout << "------------------------------------------------------" << '\n';
  std::cout << "------------------Begin Simulation--------------------" << '\n';
  std::cout << "------------------------------------------------------" << '\n';

  uint verlet_counter { 0 };      // Count successful steps since last rebinning
  long output_counter { 0 };      // Index the outputs
  bool update_neighbours { true };// Always bin neighbours on first step
  const ulong num_outputs { static_cast<ulong>(ceil(mTimeSteps/mOutFreq)) };
#ifdef PHAGE
  bool released_phage { false };
#endif
  for ( ulong tt=0; tt<=mTimeSteps; ++tt )
  {
    // check for outputting
    if ( tt%mOutFreq==0 )
    {
      updateStatus(output_counter,num_outputs);
      // std::cout << "\n=================" << '\n';
      // for ( auto &cc : mCells )
      // {
      //   std::cout << " === Cell " << cc->getID() << " ===" << '\n';
      //   for ( auto &spring : cc->getSprings() )
      //   {
      //     auto& [key,ss] = spring;
      //     std::cout << "-> " << ss.mCellB->getID()
      //               << " s: " << ss.mS << " t: " << ss.mT
      //               << " key: " << key << '\n';
      //   }
      // }

      //Output to file
      printCellsToFile(
        createFileName(
          output_counter,
          "biofilm_",
          sim_out_dir
        ),
        mCells,
        false,
        false
      );
#ifdef PHAGE
      printCellsToFile(createFileName(output_counter),mPhage,true,false);
#endif

      ++output_counter;

      if ( tt>=mTimeSteps )
      {
        break;
      }

    }

    // refresh list if every N timesteps
    if ( verlet_counter>=mGrid.mVerletUptime )
    {
      // std::cout << "Binning: verlet_counter>=mGrid.mVerletUptime" << '\n';
      update_neighbours=true;
    }
    else ++verlet_counter;

    updateOneTimeStep(update_neighbours,verlet_counter);

#ifdef PHAGE
    if ( tt*mDt*constants::baseTimeScale>0 && released_phage==false )
    {
      std::cout << "Release phage" << '\n';
      // const int num_phg { 40 };
      // for (int ii=0;ii<num_phg;++ii)
      //   mPhage.push_back(
      //     new Phage(
      //       0.5*constants::max_grid_size*gen_rand.getUniformRand(-1.0,1.0),
      //       0.5*constants::max_grid_size*gen_rand.getUniformRand(-1.0,1.0),
      //       0.1*constants::max_grid_size
      //     )
      //   );
      mPhage.push_back( new Phage(0,0,0) );
      released_phage=true;
      // for ( auto pp : mPhage )
      // {
      //   std::cout << pp->getID() << " "
      //             << 2*pp->getPos()/constants::max_grid_size << '\n';
      // }
      // exit(1);
    }
    if ( mNumberOfInfected==mCells.size() || mCells.size()==0 )
    {
      updateStatus(output_counter,num_outputs);
      printCellsToFile(
        createFileName(
          output_counter,
          "final_",
          sim_out_dir
        ),
        mCells,
        false,
        false);
      printCellsToFile(
        createFileName(
          output_counter,
          "final_",
          sim_out_dir
        ),
        mPhage,
        true,
        false);
      break;
    }
#endif
    if ( getLength()>=mTargetSize  )
    {
      printCellsToFile(
        createFileName(
          output_counter,
          "final_",
          sim_out_dir
        ),
        mCells,
        false,
        false);
        break;
    }
  }
  std::cout << '\n';
}
