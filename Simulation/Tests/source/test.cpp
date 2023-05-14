// Standard libraries
#include <iostream>
#include <cmath>
#include <vector>
#include <array>
#include <cassert>
#include <random>
#include <chrono>
#include <memory>

// Third party
#include "timer.hpp"

// Custom classes
#include "constants.hpp"         // definition of constants namespace
#include "MathUtility.hpp"
#include "particle.hpp"
#include "RodShapedBacteria.hpp"
#include "infected_class.hpp"
#include "biofilm_class.hpp"
// #include "biofilm_class.cpp"
#include "forces.hpp"
// #include "forces.cpp"
#include "chaining_susceptible_class.hpp"
// #include "chaining_susceptible_class.cpp"
// #include "test.hpp"
#include "VerletGrid.hpp"
// #include "VerletGrid.cpp"
// #include "infected_biofilm_class.hpp"

template class GridCell<RodShapedBacterium>;

void testChainingRodShapedBacteriumClass()
{
  auto chain_sus_1 = std::make_unique< ChainingRodShapedBacterium > (
    Vec3{0,0,0},0,constants::nondim_rodGrwthRtePreFac,constants::nondim_avg_div_L
  );
  std::cout << *chain_sus_1 << '\n';

  ChainingRodShapedBacterium::cellvec chained_cells;
  createDaughters(chain_sus_1,chained_cells);
  for ( auto& cell : chained_cells )
  {
    std::cout << *cell << '\n';
  }

}

void testChainingBiofilmEvolve()
{
  // auto seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
  auto seed =constants::SEED;
  std::cout << "seeding with " << seed << '\n';

  Timer t_evol_init_test;
  Biofilm<ChainingRodShapedBacterium> chaining_biofilm{};
  chaining_biofilm.evolveBiofilm(constants::conlony_size_target,seed);
  std::cout << "time: " << t_evol_init_test.elapsed() << " for "
            << chaining_biofilm.getCellsInFilm().size()
            << " cells" << '\n';
}

void testRodShapedBacteriumClass()
{

  std::cout << "Checking basic utilities for RodShapedBacterium class" << '\n';
  RodShapedBacterium cellA{};
  std::cout << "CoM at: " << cellA.mPos << '\n';
  std::cout << "Direction at: " << cellA.getOrientation() << '\n';
  std::cout << "cell area: " << cellA.getCellArea() << '\n';
  std::cout << "divide?: " << cellA.signalDivide() << '\n';
  cellA.grow(10);
  std::cout << "cell area: " << cellA.getCellArea() << '\n';
  std::cout << "divide?: " << cellA.signalDivide() << '\n';

  // test basic movement
  cellA.mForce=Vec3{1,1,0};
  cellA.mTorque=Vec3{0,0,0.5};
  cellA.mAngVel=Vec3{0,0,0.5*constants::pi};
  cellA.mVel=Vec3{1,1,0};
  cellA.move(1);
  std::cout << "CoM at: " << cellA.mPos << '\n';
  std::cout << "Direction at: " << cellA.getOrientation() << '\n';

  // test reset
  Vec3 rcm { cellA.mPos };
  Vec3 dir { cellA.getOrientation() };

  cellA.reset();
  cellA.move(0.1);
  Vec3 rcm2 { cellA.mPos };
  Vec3 dir2 { cellA.getOrientation() };
  assert ( static_cast<bool>( rcm2==rcm ) );
  assert ( static_cast<bool>( dir2==dir ) );

  std::cout << cellA << '\n';

}

void testBiofilmClassLoad()
{
  Biofilm<RodShapedBacterium> filmA{};
  std::string filename { "data/vis_biofilm_00434.txt" };
  filmA.populateCellsFromFile(filename);
  filmA.printBiofilmStateToFile( "test_output/vis_biofilm_00000.txt" );
}

void testInteractions()
{
  std::cout << "Checking basic utilities for Biofilm class" << '\n';

  /* Check Hertzian contact forces*/
  // Case 1: head to head contact
  RodShapedBacterium cellA{Vec3{0,0,0.0},0.5*constants::pi};
  GenericInfected cellC{1.6,0,0.0,0.25*constants::pi};

  RodShapedBacterium* susceptibles[]{ &cellA, &cellC };
  std::array<Vec3,2> centres{
    getVirtualContactCentres<RodShapedBacterium*>(susceptibles[0],susceptibles[1])
  };
  Vec3 cA{centres[0]};
  Vec3 cC{centres[1]};
  std::cout << "centres: cA: " << cA << " cC: " << cC << '\n';

  std::cout << "/*************Hertzian forces******************/" << '\n';
  std::array<Vec3,2> virt_centres;
  Vec3 force_B_on_A;
  getHertzianForceCellBonCellA<RodShapedBacterium*>(
    susceptibles[0],
    susceptibles[1],
    force_B_on_A,
    virt_centres
  );

  std::cout << force_B_on_A << '\n';

  std::cout << "cell A" << '\n';
  std::cout << cellA << '\n';
  std::cout << "centre: " << cA << " force: " << force_B_on_A << '\n';

  std::cout << "cell C" << '\n';
  std::cout << cellC << '\n';
  std::cout << "centre: " << cC << " force: " << -force_B_on_A << '\n';

  std::array<Vec3,2> forces;
  std::array<Vec3,2> torques;
  getHertzianForceAndTorqueBetweenCellBandCellA(
    susceptibles[0],
    susceptibles[1],
    forces,
    torques
  );

  Vec3 tau_A { torques[0] };
  Vec3 tau_C { torques[1] };
  std::cout << "tau A" << tau_A << '\n';
  std::cout << "tau C" << tau_C << '\n';

  // printVisulationPairInteraction<RodShapedBacterium>(cellA,cellC);

}

// void testBiofilmClassEvolve()
// {
//
//   // auto seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
//   auto seed { 1602366881392841974 };
//   std::cout << "seeding with " << seed << '\n';
//
//   // Create some random particles
//   std::vector< RodShapedBacterium > susceptibles;
//   std::vector< GenericInfected > infectives;
//   std::vector< Phage > phage;
//
//   // Make a bunch of randomly distributed particles
//   std::uniform_real_distribution<double> rnd_pos(-10.0, 10.0);
//
//   for ( int ii=0; ii<10; ++ii )
//   {
//     susceptibles.emplace_back(
//       rnd_pos(Particle::mGenerator),
//       rnd_pos(Particle::mGenerator),
//       0
//     );
//
//     infectives.emplace_back(
//       rnd_pos(Particle::mGenerator),
//       rnd_pos(Particle::mGenerator),
//       0
//     );
//
//     phage.emplace_back(
//       rnd_pos(Particle::mGenerator),
//       rnd_pos(Particle::mGenerator),
//       0
//     );
//
//   }
//
//   // Create a vector to control all particles
//   std::vector< RodShapedBacterium* > poly_test_vec;
//   poly_test_vec.reserve(20);
//
//   for ( auto &pp : susceptibles )
//   {
//     poly_test_vec.push_back(&pp);
//   }
//   for ( auto &pp : infectives )
//   {
//     poly_test_vec.push_back(&pp);
//   }
//
//   // Boxes must have at least a side length of the maximum distance between cells
//   GridCell<RodShapedBacterium>::box_width = 1+constants::nondim_avg_div_L;
//   std::vector< std::unique_ptr<GridCell<RodShapedBacterium>> > grid_cells;
//
//
//   std::cout << "1" << '\n';
//   allocateGrid(grid_cells,poly_test_vec);
//   std::cout << "2" << '\n';
//   populateNeighbourCells(grid_cells);        // Populate grid neighbour list
//   std::cout << "3" << '\n';
//   updateGridCellLists(grid_cells,poly_test_vec); // Bin cells to grid
//   std::cout << "4" << '\n';
//   createNeighbourLists(grid_cells,poly_test_vec);// Create cell neighbour lists
//   std::cout << "5" << '\n';
//
//   for ( auto &pp : susceptibles )
//   {
//     std::cout << "---------------------------------------------------" << '\n';
//     std::cout << pp.mId << '\n';
//     for ( auto &nn : pp.mMyNeighbourList )
//     {
//       std::cout << *nn << '\n';
//     }
//   }
//
//   Phage::phagevec all_phage;
//   RodShapedBacterium::cellvec all_susc;
//   GenericInfected::cellvec all_inf;
//   all_phage.push_back(std::make_unique<Phage>(0,0,0));
//   all_inf.push_back(std::make_unique<GenericInfected>(0,0,0));
//   all_susc.push_back(std::make_unique<RodShapedBacterium>(0,0,0));
//   InfectedBiofilm biofilm { all_susc, all_inf, all_phage };
//   biofilm.updateOneTimeStep(1e-4);
//   // std::string filename_long {
//   //   "phage_infection/data/test_biofilm.txt"
//   // };
//   // std::string filename_long {
//   //   "/media/rory/Elements/biofilm_37554/output/data/vis_biofilm_00455.txt"
//   // };
//   // std::string filename_long {
//   //   "/media/rory/Elements/biofilm_37554/output/data/vis_biofilm_00105.txt"
//   // };
//   RodShapedBacterium::mGenerator.seed(seed);
//   // biofilm::getAllForceTorque()
//   // filmC.populateCellsFromFile(filename_long);
//   // Timer t_c;
//   // filmC.getAllForceTorque();
//   // std::cout << "time: " << t_c.elapsed() << " for "
//   //           << filmC.getCellsInFilm().size()
//   //           << " cells" << '\n';
//
//   // std::cout << "The size of C is: " << filmC.size()
//   //           << " for " << filmC.getCellsInFilm().size()
//   //           << '\n';
//   // Timer t_evol;
//   // filmC.evolveBiofilm(filmC.size()+30,seed);
//   // std::cout << "time: " << t_evol.elapsed() << " for "
//   //           << filmC.getCellsInFilm().size()
//   //           << " cells" << '\n';
// }

/*
convert the parametrisation of the line from starting at the midpoint and taking
a parameter t \in [-1,1] to t_new \in [0,1] starting from one end
*/
double pourninToElbyParam(double t)
{
  return 0.5*(1+t);
}

/*
convert the parametrisation of the line from starting at the end point and taking
a parameter t \in [0,1] to t_new \in [-1,1] starting from the centre
*/
double elbyToPourninParam(double t)
{
  return 2*t -1;
}

void testVecClass()
{
  std::cout << "Checking Vec3 utilities" << '\n';

  // test initialiation methods
  // Vec3 vecB = Vec3{1,2,3};
  // Vec3 vecC = Vec3{4,5,6};
  // Vec3 vecD = Vec3{234,654,327.3};
  // assert(vecD.mY==654);

  /* test arithmetic operators */
  // test addition
  // Vec3 vecSum{vecC+vecB};
  // Vec3 check_vec{0,0,0};
  // for (int ii = 0; ii < 3; ii++) {
  //   check_vec[ii] = rndvec[ii] + rndvec2[ii];
  //   assert(check_vec[ii] == vecSum[ii]);
  // }
  //
  // assert(check_vec == vecSum);

  // test unary operator-
  // for (int ii = 0; ii < 3; ii++) {
  //   check_vec[ii] = -(rndvec[ii] + rndvec2[ii]);
  //   assert(check_vec[ii] == (-vecSum)[ii]);
  // }

  // test binary minus
  // for(int ii=0; ii< 3;++ii)
  // {
  //   assert(vecSum[ii]-vecC[ii] == vecB[ii]);
  // }

  // test +=
  // std::cout << vecB << '\n';
  // vecB += vecC;
  // assert( vecB == vecSum );
  // std::cout << vecB << '\n';
  // vecSum -= vecC;
  // assert ( static_cast<bool>(vecSum == Vec3{1,2,3}) );

  /* Test vector operators */
  // double resDot1{vecB.dot(vecC)};
  // std::cout << "vecB.vecC " << resDot1 << '\n';

  Vec3 vecB, vecC;
  Vec3 cA{1,2,3}, cB{3,2,1};
  assert( eulerianDot(cA,cB) == 10 );
  assert( cA.norm() == cB.norm() );
  assert( cA.norm() == sqrt(14) );

  Vec3 cr_vec{vecB.cross(vecC)};
  // std::cout << "vecB x vecC " << cr_vec << '\n';

  Vec3 vec_i{1,0,0}; Vec3 vec_j{0,1,0}; Vec3 vec_k{0,0,1};
  std::cout << vec_i.cross(vec_j) << " " << eulerianCross(vec_i,vec_j) << '\n';
  assert(vec_i.cross(vec_j) ==  vec_k );
  assert(vec_j.cross(vec_i) == -vec_k );
  assert(vec_k.cross(vec_i) ==  vec_j );
  // std::cout << vec_i.cross(vec_k) << '\n';

  assert(angularSep(vec_i,vec_k) / constants::pi == 0.5);
  // std::cout << "angle between i and k is "
  //           << angularSep(vec_i,vec_k) / constants::pi
  //           << " pi radians" << '\n';

  std::cout << "Check parallel closest distance" << '\n';
  Vec3 x1{0,0,0};
  Vec3 a1{1,0,0};
  Vec3 x2{4,0,0};
  Vec3 a2{1,0,0};
  double s_star,t_star;

  // Check case of no overlap - only extrema should be returned
  closestApproachLineSegmentsPournin(x1,a1,x2,a2,s_star,t_star);
  assert(s_star==1); assert(t_star==-1);
  closestApproachLineSegmentsPournin(x1,a1,x2,-a2,s_star,t_star);
  assert(s_star==1); assert(t_star==1);
  closestApproachLineSegmentsPournin(x1,-a1,x2,a2,s_star,t_star);
  assert(s_star==-1); assert(t_star==-1);
  closestApproachLineSegmentsPournin(x1,-a1,x2,-a2,s_star,t_star);
  assert(s_star==-1); assert(t_star==1);

  // Check the case of partial overlap - centres are still outside the other segment
  a1 = Vec3{3.1,0,0};
  closestApproachLineSegmentsPournin(x1,a1,x2,a2,s_star,t_star);
  assert(  fabs( s_star - (1-0.05/3.1) ) < 1e-15  );
  assert(  fabs( t_star + 0.95 ) < 1e-15 );

  // Check the case of centres projecting into the other segment
  a1 = Vec3{5,0,0}; a2 = Vec3{1,0,0}; // segment 2 inside segment 1
  closestApproachLineSegmentsPournin(x1,a1,x2,a2,s_star,t_star);
  assert(s_star==0.8); assert(t_star==0);

  a1 = Vec3{1,0,0}; a2 = Vec3{10,0,0}; // segment 1 inside segment 2
  closestApproachLineSegmentsPournin(x1,a1,x2,a2,s_star,t_star);
  assert(s_star==0); assert(t_star==-0.4);
  closestApproachLineSegmentsPournin(x1,a1,x2,-a2,s_star,t_star);
  assert(s_star==0); assert(t_star==0.4);

  // test from [2]
  Vec3 A{0,0,0}; Vec3 B{1,2,1}; Vec3 C{1,0,0}; Vec3 D{2,1,0};
  x1 = 0.5 * (A + B); x2 = 0.5 * (C + D); a1 = 0.5 * (B-A); a2 = 0.5 * (D-C);
  closestApproachLineSegmentsPournin(x1,a1,x2,a2,s_star,t_star);
  Vec3 p{x1+s_star*a1};
  Vec3 q{x2+t_star*a2};
  assert( ( eulerianDot(p-q,p-q) - 0.8333 ) < 1e-3 );

  // test from https://www.geometrictools.com/Documentation/DistanceLine3Line3.pdf
  A = Vec3{-1.0896,  9.7235e-7, 0.0};
  B = Vec3{ 0.9122, -9.4370e-7, 0.0 };
  C = Vec3{-0.9001,  9.0671e-7, 0.0 };
  D = Vec3{ 1.0731, -9.8186e-7, 0.0};
  x1 = 0.5 * (A + B); x2 = 0.5 * (C + D); a1 = 0.5 * (B-A); a2 = 0.5 * (D-C);
  closestApproachLineSegmentsPournin(x1,a1,x2,a2,s_star,t_star);
  p = x1+s_star*a1;
  q = x2+t_star*a2;
  assert( (eulerianNorm(p-q) - 1.1575e-7) < 1e-10 );

  // Sunday test case
  A = Vec3{ 0.779990, 0.611925, -0.227031 };
  B = Vec3{ 0.532153, 0.857246, -0.101024 };
  C = Vec3{-0.212773, 0.350915, -0.495572 };
  D = Vec3{ 0.118815, 0.022495, -0.664266 };
  x1 = 0.5 * (A + B); x2 = 0.5 * (C + D); a1 = 0.5 * (B-A); a2 = 0.5 * (D-C);
  closestApproachLineSegmentsPournin(x1,a1,x2,a2,s_star,t_star);
  p = x1+s_star*a1;
  q = x2+t_star*a2;
  std::cout << "eulerianNorm(p-q): " << eulerianNorm(p-q)
            <<  " s_star: " << s_star << " t_star " << t_star << '\n';
  assert( (eulerianNorm(p-q) - 0.982924) < 1e-3 );

  // Check intersection
  double delta{ 0.25*1e-4 };
  double eps{ sqrt(delta) };
  double phi{ 1e-5 };
  A = Vec3{0,0,0}; B = Vec3{ 1,0,0 };
  C = Vec3{ -eps, phi+delta, 0 }; D = Vec3{ eps, phi-delta, 0 };
  x1 = 0.5 * (A + B); x2 = 0.5 * (C + D); a1 = 0.5 * (B-A); a2 = 0.5 * (D-C);
  closestApproachLineSegmentsPournin(x1,a1,x2,a2,s_star,t_star);
  p = x1+s_star*a1;
  q = x2+t_star*a2;
  assert( ( pourninToElbyParam(s_star) - 0.002 ) < 1e-15 );
  assert( ( pourninToElbyParam(t_star) - 0.700 ) < 1e-15 );
  std::cout <<  "eulerianNorm(p-q): " << eulerianNorm(p-q)
            <<  " s_star: " << s_star
            <<  " t_star "  << t_star
            << '\n';

  // Check code test case
  // double d_0{1};
  double L{2.00433};
  x1 = Vec3{1.09325,1.10222,0};
  a1 = 0.5*L*Vec3{0.948984,0.315325,0};
  x2 = Vec3{1.68993,0.537289,0};
  a2 = 0.5*L*Vec3{0.879442,0.476006,0};
  std::cout << "a1.a2: " << eulerianDot(a1,a2) << '\n';
  closestApproachLineSegmentsPournin(x1,a1,x2,a2,s_star,t_star);
  p = x1+s_star*a1;
  q = x2+t_star*a2;
  assert( (eulerianNorm(p-q) - 0.606058) < 1e-15 );


  std::cout << "Test of vec class completed successfully" << '\n';
}

int main(int argc, char const *argv[])
{
  std::cout << "gamma: " << constants::nondim_rodGrwthRtePreFac << "\n"
            << "dt: "    << constants::nondim_dt << '\n';
  testVecClass();
  testRodShapedBacteriumClass();
  testInteractions();
  // testBiofilmClassEvolve();
  // testBiofilmClassLoad();
  testChainingRodShapedBacteriumClass();
  testChainingBiofilmEvolve();
  return 0;
}
