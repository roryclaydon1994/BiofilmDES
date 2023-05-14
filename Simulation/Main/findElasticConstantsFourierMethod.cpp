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
#include <omp.h>
#include <filesystem>
#include <functional>
// //#include <execution>

// Custom classes
#include "constants.hpp"         // definition of constants namespace
#include "MathUtility.hpp"
#include "RandUtil.hpp"
#include "IBacterium.hpp"
#include "SphericalBacteria.hpp"
#include "RodShapedBacteria.hpp"
#include "VerletGrid.hpp"
#include "forces.hpp"
#include "IO.hpp"
#include "PolyBiofilm.hpp"

constexpr int nlines { 5 };
enum col_names {
  straight_channel_initial_hertzian,
  straight_channel_initial_bend,
  straight_channel_initial_ext,
  annulus_channel_initial_hertzian,
  annulus_channel_initial_bend,
  annulus_channel_initial_ext,
  straight_channel_relaxed_hertzian,
  straight_channel_relaxed_bend,
  straight_channel_relaxed_ext,
  annulus_channel_relaxed_hertzian,
  annulus_channel_relaxed_bend,
  annulus_channel_relaxed_ext,
  normalised_initial_delta_E,
  normalised_relaxed_delta_E
};

#ifdef CHAINING

using BVec = std::vector<IBacterium*>;
using Earray = std::array<double,14>;

void getDeltaEnergy(
  Earray &energies,
  uint delta_idx,
  uint sh, uint sb, uint se,
  uint ah, uint ab, uint ae,
  double area
)
{
  double straight_energy { energies[sh]+energies[sb]+energies[se] };
  double annulus_energy  { energies[ah]+energies[ab]+energies[ae] };
  energies[delta_idx]=(annulus_energy-straight_energy)/area;
}

double getOverlap(
  IBacterium* A,
  IBacterium* B
)
{
  double s1,t1,s2,t2; // Line parameter values at which they are closest
  Vec3 c1,c2;         // centres of the closet approach on each line resp.

  // COMs
  Vec3 pos_A { A->getPos() };
  Vec3 pos_B { B->getPos() };

  // find the direction Segment of S1 and S2
  Vec3 v1 = 0.5*A->getLength()*A->getOrientation();
  Vec3 v2 = 0.5*B->getLength()*B->getOrientation();

  // Check line was parallel
  closestApproachLineSegmentsParallel(pos_A,v1,pos_B,v2,s1,t1,s2,t2);
  // if ( par )
  // {
  //   std::cout << "Cells should not be parallel here" << '\n';
  //   // exit(89);
  // }
  c1 = pos_A+v1*s1; // centre of the virtual sphere on cell A
  c2 = pos_B+v2*t1; // centre of the virtual sphere on cell B

  Vec3 cv = c1-c2;
  double sep=cv.norm();
  double sigma { A->getRadius() + B->getRadius() };
  double overlap { std::max(sigma-sep,0.0) }; // Cell overlap
  return overlap;
}

double getMaxPairOverlap(IBacterium* A, IBacterium* B)
{

#ifdef ADHESION
  std::cout << "ADHESION not supported at present" << '\n';
  exit(12);
#endif

  // double energy {0.0};// Energy due to an overlap
  double s1,t1,s2,t2; // Line parameter values at which they are closest
  Vec3 c1,c2;         // centres of the closet approach on each line resp.

  // COMs
  Vec3 pos_A { A->getPos() };
  Vec3 pos_B { B->getPos() };

  // find the direction Segment of S1 and S2
  Vec3 v1 = 0.5*A->getLength()*A->getOrientation();
  Vec3 v2 = 0.5*B->getLength()*B->getOrientation();

  // Check line was parallel
  bool par = closestApproachLineSegmentsParallel(pos_A,v1,pos_B,v2,s1,t1,s2,t2);

  c1 = pos_A+v1*s1; // centre of the virtual sphere on cell A
  c2 = pos_B+v2*t1; // centre of the virtual sphere on cell B
  double overlap { std::max(1-(c1-c2).norm(),0.0) };

  // Parallel segements so two points were returned
  if ( par )
  {
    c1 = pos_A+v1*s2; // centre of the virtual sphere on cell A
    c2 = pos_B+v2*t2; // centre of the virtual sphere on cell B
    double tm_ov { std::max(1-(c1-c2).norm(),0.0) };
    overlap=std::max(overlap,tm_ov);
  }
  return overlap;
}

void createLogFile(
  std::string run_dir,
  uint ntrials,
  double radius,
  double width
)
{
  std::filesystem::create_directories(sim_out_dir);
  std::string filename=run_dir+'/'+"LcDist.log";
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
  headers.push_back("RodAspectRatio");
  values.push_back(RodShapedBacterium::mAvgDivLen/constants::rodSpheroDiam);

  headers.push_back("RodModE");
  values.push_back(RodShapedBacterium::mRodModE);

  headers.push_back("RodGrowthRate");
  values.push_back(RodShapedBacterium::mAvgGrwthRate);

  headers.push_back("Kappa");
  values.push_back(RodShapedBacterium::mKappa);

  headers.push_back("BendRig");
  values.push_back(RodShapedBacterium::mBendRig);

  headers.push_back("LinkingProb");
  values.push_back(RodShapedBacterium::mLinkingProb);

  headers.push_back("ForceThresh");
  values.push_back(RodShapedBacterium::mForceThresh);

  headers.push_back("NTrials");
  values.push_back(ntrials);

  headers.push_back("annulus_radius");
  values.push_back(radius);

  headers.push_back("annulus_width");
  values.push_back(width);

  for ( uint ii=0; ii<headers.size(); ++ii )
  {
    if ( ii<headers.size()-1 )
      out_file << headers[ii] << '\t';
    else
      out_file << headers[ii] << '\n';
  }

  for ( uint ii=0; ii<values.size(); ++ii )
  {
    if ( ii<headers.size()-1 )
      out_file << values[ii] << '\t';
    else
      out_file << values[ii] << '\n';
  }
  out_file.close();
  std::cout << "\nSaved log file to " << filename << '\n';
}

void saveDists(
  std::string run_dir,
  const std::vector<std::array<double,14>> &energies
)
{
  std::string filename=run_dir+'/'+"energies.dat";
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

  out_file << "straight_channel_initial_hertzian" << '\t';
  out_file << "straight_channel_initial_bend"     << '\t';
  out_file << "straight_channel_initial_ext"      << '\t';
  out_file << "annulus_channel_initial_hertzian"  << '\t';
  out_file << "annulus_channel_initial_bend"      << '\t';
  out_file << "annulus_channel_initial_ext"       << '\t';
  out_file << "straight_channel_relaxed_hertzian" << '\t';
  out_file << "straight_channel_relaxed_bend"     << '\t';
  out_file << "straight_channel_relaxed_ext"      << '\t';
  out_file << "annulus_channel_relaxed_hertzian"  << '\t';
  out_file << "annulus_channel_relaxed_bend"      << '\t';
  out_file << "annulus_channel_relaxed_ext"       << '\t';
  out_file << "normalised_initial_delta_E"        << '\t';
  out_file << "normalised_relaxed_delta_E"        << '\n';

  for ( uint ii=0; ii<energies.size(); ++ii )
  {
    for ( uint jj=0; jj<energies[ii].size(); ++jj )
      if (jj==energies[ii].size()-1)
        out_file << energies[ii][jj] << '\n';
      else
        out_file << energies[ii][jj] << '\t';
  }
  out_file.close();

  std::cout << "Saved log file to " << filename << '\n';
}

template <class T>
void addLinksToLine(T &cells)
{
  for ( uint ii=1; ii<cells.size(); ++ii )
  {
    double link_prob { gen_rand.getUniformRand(0,1) };
    // std::cout << "link_prob: " << link_prob << '\n';
    if ( link_prob <= RodShapedBacterium::mLinkingProb )
    {
      cells[ii-1]->setUpperLink(cells[ii]);
      cells[ii]->setLowerLink(cells[ii-1]);
    }
  }
}

double getMaximumOvlerap(BVec &cells)
{
  double overlap{0.0};
  #pragma omp parallel for shared(cells) \
          schedule(static) default(none) \
          reduction(max:overlap)
  for ( uint ii=0;ii<cells.size();++ii )
  {
    IBacterium* cell { cells[ii] };
    for ( auto &neighbour_cell : cell->getNeighbourList() )
    {
      const bool different_cells {
        ( cell->getID()!=neighbour_cell->getID() )
        ||
        ( cell->getMyType()!=neighbour_cell->getMyType() )
      };
      if ( different_cells )
      {
        double tmp_overlap { getMaxPairOverlap(cell,neighbour_cell) };
        // printf("%f\n",tmp_overlap);
        overlap=std::max(overlap,tmp_overlap);
      }
    }
  }
  std::cout << "overlap: " << overlap << '\n';
  return overlap;// double counting
}

double getHertzianEnergy(BVec &cells)
{
  double energy{0.0};
  #pragma omp parallel for shared(cells) \
          schedule(static) default(none) \
          reduction(+:energy)
  for ( uint ii=0;ii<cells.size();++ii )
  {
    IBacterium* cell { cells[ii] };
    for ( auto &neighbour_cell : cell->getNeighbourList() )
    {
      const bool different_cells {
        ( cell->getID()!=neighbour_cell->getID() )
        ||
        ( cell->getMyType()!=neighbour_cell->getMyType() )
      };
      if ( different_cells )
      {
        // double ov { getMaxPairOverlap(cell,neighbour_cell) };
        // printf("%f for %f\n",tmp_energy,pow(ov,2.5));
        // double tmp_energy { 0.4*ov*ov*sqrt(ov) };
        double tmp_energy { getPairHertzianEnergy(cell,neighbour_cell) };
        energy+=tmp_energy;
      }
    }
  }
  // std::cout << "energy: " << energy << '\n';
  return 0.5*energy;// double counting
}

void getSpringEnergy(BVec &cells,double &bend_energy, double &ext_energy)
{
  bend_energy=0.0; ext_energy=0.0;
  for ( auto cell : cells )
  {
    if ( cell->getUpperLink() )
    {
      double tmp_bend_energy,tmp_ext_energy;
      getSpringEnergy(
        cell,cell->getUpperLink(),
        tmp_bend_energy,tmp_ext_energy
      );
      bend_energy+=tmp_bend_energy;
      ext_energy+=tmp_ext_energy;
    }
  }
}

double getTotalEnergy(BVec &cells)
{
  double bend_energy,ext_energy;
  getSpringEnergy(cells,bend_energy,ext_energy);
  double steric_energy { getHertzianEnergy(cells) };
  return steric_energy+bend_energy+ext_energy;
}

void logEnergies(PolyBiofilm &pb, Earray &energies, int h, int b, int e)
{
  double bend_energy,ext_energy;
  pb.mGrid.updateVerletList(pb.mCells);
  energies[h]=getHertzianEnergy(pb.mCells);
  getSpringEnergy(pb.mCells,bend_energy,ext_energy);
  energies[b]=bend_energy;
  energies[e]=ext_energy;
}

void createStraight(
  double channel_length,
  double height,
  BVec &cells
)
{
  double l_0 {
    0.5*( RodShapedBacterium::mAvgDivLen-2*RodShapedBacterium::mRadius )
  };
  std::vector<RodShapedBacterium*> cells_in_line;

  double length {
    gen_rand.getUniformRand(l_0,RodShapedBacterium::mAvgDivLen )
  };

  double start_x { -0.5*channel_length+0.5+0.5*length };
  RodShapedBacterium* rod = new RodShapedBacterium{
    start_x,height,0,    // position
    0,constants::pi*0.5, // angles
    0.0,                 // growth rate
    length               // random length (uniform for now)
  };
  cells_in_line.push_back(rod);
  double total_length { length+1 };
  while( isclose(total_length,channel_length)==false )
  {
    double space_left { channel_length-total_length };
    double new_length {
      gen_rand.getUniformRand(l_0,RodShapedBacterium::mAvgDivLen )
    };
    new_length=std::min(new_length,space_left-1);
    if ( space_left>=(l_0+1) && space_left<=(RodShapedBacterium::mAvgDivLen+1) )
    {
      new_length=space_left-1;
    }
    else if (
      (RodShapedBacterium::mAvgDivLen+1)<space_left && space_left<2*(RodShapedBacterium::mAvgDivLen+1)
    )
    {
      new_length=std::min(RodShapedBacterium::mAvgDivLen,space_left-2-l_0);
    }
    else if ( space_left<(l_0+1) )
    {
      std::cout << "last x " << cells_in_line.back()->mPos.x+0.5*cells_in_line.back()->mLength+0.5 << '\n';
      std::cout << "fail " << space_left << '\n';
      exit(23);
    }

    double x {
      cells_in_line.back()->mPos.x+0.5*cells_in_line.back()->mLength+
      0.5*new_length+1
    };
    RodShapedBacterium* new_rod = new RodShapedBacterium{
      x,height,0,          // position
      0,constants::pi*0.5, // angles
      0.0,                 // growth rate
      new_length           // random length (uniform for now)
    };

    cells_in_line.push_back(new_rod);
    total_length+=new_length+1;
  }


  double final_length {
    std::accumulate(
      cells_in_line.begin(),cells_in_line.end(),
      0.0,
      [] (double lhs, auto cb) -> double { return lhs+cb->getLength()+1; }
    )
  };

  if ( final_length<channel_length-1e3 )
  {
    std::cout << "Error! Final length is too short!" << '\n';
  }
  addLinksToLine(cells_in_line);
  cells.insert(cells.end(),cells_in_line.begin(),cells_in_line.end());
}

void setStraightCondition(
  const double radius,
  const double width,
  BVec &cells
)
{
  // Set the centre to centre spacing of the outer cells s.t. width is correct
  double spacings { (width-1.0)/(2.0*nlines) };
  for ( int ii=-nlines; ii<=nlines; ++ii )
  {
    BVec cells_in_line;
    double height { ii*spacings };
    createStraight(radius,height,cells_in_line);
    cells.insert(cells.end(),cells_in_line.begin(),cells_in_line.end());
  }
}

void relaxSystem(
  double length,
  double width,
  BVec &interior_cells,
  BVec &boundary_cells
)
{
  std::cout << "=======================================================" << '\n';
  // std::unordered_map<uint,bool> onBoundary;
  // onBoundary.reserve(pb.mCells.size());
  // for ( auto &cell : interior_cells ) onBoundary[cell->getID()]=false;
  // for ( auto &cell : boundary_cells ) onBoundary[cell->getID()]=true;
  // for ( auto it=interior_cells.begin(); it!=interior_cells.end(); )
  // {
  //   if ( onBoundary[(*it)->getID()] ) it=interior_cells.erase(it);
  //   else ++it;
  // }
  // std::cout << "bcs: " << boundary_cells.size() << '\n';

  // pb.mOutFreq=static_cast<uint>(
  //    ceil( ( 0.02/constants::baseTimeScale )/pb.mDt)
  // );
  // pb.mTimeSteps=static_cast<ulong>(
  //   ceil( ( 0.1/constants::baseTimeScale )/pb.mDt)
  // );
  // pb.runSim();
  BVec all_cells;
  all_cells.insert(all_cells.end(),interior_cells.begin(),interior_cells.end());
  all_cells.insert(all_cells.end(),boundary_cells.begin(),boundary_cells.end());
  GridCells<IBacterium> grid { all_cells,constants::box_width };
  grid.updateVerletList(all_cells);
  double old_energy { getTotalEnergy(interior_cells) };
  // double original_energy { old_energy };
  double delta_energy;
  uint counter { 0 };
  do
  {
    if ( counter>1e4 )
    {
      std::cout << "Error, too many steps to converge." << '\n';
      exit(34);
    }
    polyInteractParticles(all_cells);
    // applyBoundaryCondition(boundary_cells,length,width);
    for ( auto &cell : interior_cells )
    {
      cell->move(constants::nondim_dt);
      double moved {
        dot2( cell->getLoggedPos()-cell->getPos() )
      };
      if ( moved>dot2( grid.mVerletMoveR ) )
      {
        // std::cout << "Binning: moved too far" << '\n';
        grid.updateVerletList(all_cells);
      }
      cell->reset();
    }
    double new_energy { getTotalEnergy(interior_cells) };
    delta_energy=fabs( new_energy-old_energy );
    if ( counter%1000==0 )
    {
      std::cout << " N: " << new_energy
                << " O: " << old_energy
                << " D: " << delta_energy
                << " D>e: " << (delta_energy>1e-3)
                << '\n';
      getMaximumOvlerap(interior_cells);
    }
    old_energy=new_energy;
    ++counter;
  } while ( delta_energy>1e-3 );
  std::cout << "=======================================================" << '\n';
}

// void relaxAnnulus(PolyBiofilm &pb, double length, double width)
// {
//   std::cout << "Annulus" << '\n';
//   double old_energy { getHertzianEnergy(pb.mCells) };
//   double delta_energy { 1e10 };
//   uint counter { 0 };
//   while ( delta_energy>1e-4 )
//   {
//     polyInteractParticles(pb.mCells);
//     for ( auto &cell : pb.mCells ) cell->move(pb.mDt);
//     double new_energy { getHertzianEnergy(pb.mCells) };
//     delta_energy=fabs( (new_energy-old_energy) / old_energy );
//     std::cout << " N: " << new_energy
//               << " O: " << old_energy
//               << " D: " << delta_energy
//               << '\n';
//     old_energy=new_energy;
//     ++counter;
//     if ( counter==10 )
//     {
//       std::cout << "update cells" << '\n';
//       pb.mGrid.updateVerletList(pb.mCells);
//       counter=0;
//     }
//   }
// }

void runOrientationalStiffnessSim(
  const double radius,
  const double width,
  std::vector<std::array<double,14>> &energies,
  const uint trial_num,
  const std::string save_dir
)
{
  constexpr bool relax { true };

  std::stringstream ss;
  ss << std::setw(5) << std::setfill('0') << trial_num;

  // Hold a list of cells along the boundary. These cells will not move.
  BVec boundary_cells;
  setStraightBoundaryCondition(
    2*constants::pi*radius,width,
    boundary_cells
  );
  if ( trial_num==0 )
    printCellsToFile(
      createFileName(
        "straight_channel_boundary_"+ss.str(),
        "biofilm_",
        save_dir
      ),
      boundary_cells,
      false,
      false
    );

  // Set the initial conditions for the straight channel
  BVec straight_channel_cells;
  setStraightCondition(
    2*constants::pi*radius,width,
    straight_channel_cells
  );
  if ( trial_num==0 )
    printCellsToFile(
      createFileName(
        "straight_channel_initial_"+ss.str(),
        "biofilm_",
        save_dir
      ),
      straight_channel_cells,
      false,
      false);

  // Set up biofilm
  PolyBiofilm pb_straight {
    straight_channel_cells
  };
  logEnergies(
    pb_straight,energies[trial_num],straight_channel_initial_hertzian,
    straight_channel_initial_bend,straight_channel_initial_ext
  );

  std::cout << "Initial straight max overlap" << '\n';
  pb_straight.mGrid.updateVerletList(pb_straight.mCells);
  getMaximumOvlerap(pb_straight.mCells);
  // std::cout << "turned off relaxation" << '\n';
  if ( relax )
    relaxSystem(
      2*constants::pi*radius,
      width,
      straight_channel_cells,
      boundary_cells
    );
  // for ( auto cell : straight_channel_cells )
  //   cell->printToFile(std::cout);
  std::cout << "Relaxed straight max overlap" << '\n';
  pb_straight.mGrid.updateVerletList(pb_straight.mCells);
  getMaximumOvlerap(pb_straight.mCells);

  if ( trial_num==0 )
    printCellsToFile(
      createFileName(
        "straight_channel_relaxed_"+ss.str(),
        "biofilm_",
        save_dir
      ),
      straight_channel_cells,
      false,
      false);

  // std::cout << "Second energy" << '\n';
  logEnergies(
    pb_straight,energies[trial_num],straight_channel_relaxed_hertzian,
    straight_channel_relaxed_bend,straight_channel_relaxed_ext
  );

  // Set the centre to centre spacing of the outer cells s.t. width is correct
  double spacings { (width-1.0)/(2.0*nlines) };
  // std::cout << "spacings: " << spacings << '\n';

  // Set the initial conditions
  // Create a lambda function with R_upper and R_lower
  boundary_cells.clear();
  double R_upper { radius+spacings*nlines+1 };
  double R_lower { radius-spacings*nlines-1 };
  setAnnulusConditionBoundaryCondition(
    {R_upper,R_lower},
    boundary_cells
  );
  if ( trial_num==0 )
    printCellsToFile(
      createFileName(
        "annulus_channel_boundary_"+ss.str(),
        "biofilm_",
        save_dir
      ),
      boundary_cells,
      false,
      false
    );


  BVec annulus_cells;
  setAnnulusCondition(
    radius,
    width,
    annulus_cells
  );

  // Set up biofilm
  PolyBiofilm pb_annulus {
    annulus_cells
  };
  logEnergies(
    pb_annulus,energies[trial_num],annulus_channel_initial_hertzian,
    annulus_channel_initial_bend,annulus_channel_initial_ext
  );
  if ( trial_num==0 )
    printCellsToFile(
      createFileName(
        "annulus_channel_initial_"+ss.str(),
        "biofilm_",
        save_dir
      ),
      pb_annulus.mCells,
      false,
      false);

  std::cout << "Initial annulus max overlap" << '\n';
  pb_annulus.mGrid.updateVerletList(pb_annulus.mCells);
  getMaximumOvlerap(pb_annulus.mCells);
  // std::cout << "turned off relaxation" << '\n';
  if ( relax )
    relaxSystem(
      radius,
      width,
      annulus_cells,
      boundary_cells
    );
  std::cout << "Relaxed annulus max overlap" << '\n';
  pb_annulus.mGrid.updateVerletList(pb_annulus.mCells);
  getMaximumOvlerap(pb_annulus.mCells);

  if ( trial_num==0 )
    printCellsToFile(
      createFileName(
        "annulus_channel_relaxed_"+ss.str(),
        "biofilm_",
        save_dir
      ),
      pb_annulus.mCells,
      false,
      false);

  logEnergies(
    pb_annulus,energies[trial_num],annulus_channel_relaxed_hertzian,
    annulus_channel_relaxed_bend,annulus_channel_relaxed_ext
  );

  double area { 2*constants::pi*width*radius };
  getDeltaEnergy(
    energies[trial_num],
    normalised_initial_delta_E,
    straight_channel_initial_hertzian,
    straight_channel_initial_bend,
    straight_channel_initial_ext,
    annulus_channel_initial_hertzian,
    annulus_channel_initial_bend,
    annulus_channel_initial_ext,
    area
  );
  getDeltaEnergy(
    energies[trial_num],
    normalised_relaxed_delta_E,
    straight_channel_relaxed_hertzian,
    straight_channel_relaxed_bend,
    straight_channel_relaxed_ext,
    annulus_channel_relaxed_hertzian,
    annulus_channel_relaxed_bend,
    annulus_channel_relaxed_ext,
    area
  );
}

void findOrientationalScaling(
  double radius,
  double width,
  double linking_prob,
  double kappa,
  double bend_rig
)
{
  constexpr double force_thresh { 1e6 };
  RodShapedBacterium::mAvgGrwthRate = 0.0;

  initialiseChainingParameters(
    kappa,
    bend_rig,
    linking_prob,
    force_thresh
  );

  // Number of trials per radius of curvature
  constexpr uint n_trials_per_curv { 3 };

  std::stringstream ss;
  ss << "radius_" << radius
     << "_width_" << width
     << "_linking_prob_" << linking_prob
     << "_kappa_" << kappa
     << "_B_" << bend_rig;
  std::string run_dir { sim_out_dir + "/" + ss.str() + '/' };
  std::filesystem::create_directories(run_dir);

  std::vector<std::array<double,14>> energies(n_trials_per_curv);

  for ( uint jj=0; jj<n_trials_per_curv; ++jj )
  {
    gen_rand.setRandomSeed();
    runOrientationalStiffnessSim(
      radius,
      width,
      energies,
      jj,
      run_dir
    );
  }

  createLogFile(
    run_dir,
    n_trials_per_curv,
    radius,
    width
  );
  saveDists(
    run_dir,
    energies
  );
}

#endif

int main(int argc, char const *argv[])
{
#ifndef CHAINING
  std::cout << "Please recompile w/ CHAINING" << '\n';
  exit(2);
#else
  if ( argc==7 )
  {
    const std::string run_dir {            argv[1]  }; // Run directory
    const double radius       { std::stod( argv[2]) }; // Radius of curvature
    const double width        { std::stod( argv[3]) }; // Channel width
    const double linking_prob { std::stod( argv[4]) }; // Linking probability
    const double kappa        { std::stod( argv[5]) }; // Spring constant
    const double bend_rig     { std::stod( argv[6]) }; // Bending rigidity

    sim_out_dir += "/" + run_dir + "/";
    findOrientationalScaling(radius,width,linking_prob,kappa,bend_rig);
  }
  else
  {
    std::cout << "Expected 5 command line arguments! Received " << argc-1 << '\n';
    std::cout << "Example usage:\n"
              << "./main.out run_dir radius width linking_prob kappa bend_rig " << '\n';
    exit(EXIT_FAILURE);
  }
  return 0;
  #endif
}
