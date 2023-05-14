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
// //#include <execution>

// Custom classes
#include "constants.hpp"         // definition of constants namespace
#include "MathUtility.hpp"
#include "particle.hpp"
#include "RodShapedBacteria.hpp"
#include "chaining_susceptible_class.hpp"
#include "infected_class.hpp"
#include "Phage.hpp"
#include "VerletGrid.hpp"
#include "forces.hpp"
#include "HyperParams.hpp"
#include "IO.hpp"
#include "biofilm_class.hpp"

std::vector< RodShapedBacterium >     sus;
std::vector< GenericInfected > inf;
std::vector< Phage >    phg;
IndexGrid tg;
HyperParams hyperParams;

void setUp()
{
  // Read in from file
  populateCellsFromFile(
    "/home/rory/PhD_sync/BiofilmPhageDES/GeneratedOutput/InitialConditions/init_large_1.txt",
    sus,
    true
  );
  // Set the grid and update the HyperParams
  tg.createCLL(sus,inf,phg);
  return;
}

void testSusBiofilm()
{
  setUp();
  Biofilm<RodShapedBacterium> biofilm1{};
  Biofilm<RodShapedBacterium> biofilm2{sus};
  biofilm2.updateOneTimeStep(0.01);
  biofilm1.evolveBiofilm(100);
}

void testChainingBiofilm()
{
  Biofilm<ChainingRodShapedBacterium> biofilm{};
  std::cout << " ----- Check Evolve ------ " << '\n';
  biofilm.evolveBiofilm(50);
}

void testChainingDivision()
{
  std::vector<ChainingRodShapedBacterium> ss;
  ss.emplace_back
  (
    0, 0, 0,
    0, constants::pi*0.5,
    constants::nondim_rodGrwthRtePreFac,
    constants::nondim_avg_div_L,
    BIG, BIG
  );
  Biofilm<ChainingRodShapedBacterium> biofilm{ss};
  for ( auto &cell : biofilm.mCellList )
  {
    std::cout << cell << '\n';
  }
  std::cout << " ----- Check divisions ------ " << '\n';
  biofilm.handleAllDivisionEvents();
}

int main(int argc, char const *argv[])
{
  testChainingDivision();
  testChainingBiofilm();
  // testSusBiofilm();
  return 0;
}
