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
#include "IndexGrid.hpp"
#include "HyperParams.hpp"

HyperParams hyperParams;

std::mt19937 generator {
  static_cast<unsigned long>(
    std::chrono::high_resolution_clock::now().time_since_epoch().count()
  )
};

// std::vector< Particle* >        par(1);
constexpr long NUM { 10 };
IndexGrid tg;
std::vector< RodShapedBacterium >       sus(NUM);
std::vector< GenericInfected > inf(NUM);
std::vector< Phage >    phg(NUM);


std::vector< RodShapedBacterium > sus_check;

inline bool closeParticles( RodShapedBacterium &A, RodShapedBacterium &B )
{
  return dot2(B.mPos-A.mPos)<pow(1+0.5*( A.mLength+B.mLength ),2);
}

template <class T>
bool checkTgInteractions(
  uint my_id,
  uint jj,
  std::vector<T> &cell_begin,
  std::vector<T> &cell_end
)
{
  uint index_start = cell_begin[my_id];
  if ( index_start != BIG )
  {
    uint index_end = cell_end[my_id];
    for ( uint np_id=index_start; np_id<index_end; ++np_id )
    {
      if ( my_id!=np_id && jj==np_id && closeParticles(sus_check[my_id],sus_check[np_id]) )
      {
        return true;
      }
    }
  }
  return false;
}

// void neighbourCheck()
// {
//   const long N { sus_check.size() };
//   assert( N>0 );
//   for ( long ii=0; ii<N; ++ii )
//   {
//     for ( long jj=0; jj<N; ++jj )
//     {
//       bool interaction_flag {
//         checkTgInteractions(
//           ii,jj,
//           tg.mCellBegin[IndexGrid::Types::SUS],
//           tg.mCellEnd[IndexGrid::Types::SUS],
//         )
//       };
//       if (ii==jj) continue;
//       if ( dot2(sus_check[jj].mPos-sus_check[ii].mPos)<pow(1+0.5*( sus_check[ii].mLength+sus_check[jj].mLength ),2) )
//       {
//         // ii_gp = getCellIndex(sus_check[ii].mPos);
//         // jj_gp = getCellIndex(sus_check[jj].mPos);
//         // std::cout << ii_gp << " " << jj_gp << '\n';
//       }
//     }
//   }
// }

void setUp()
{
  std::normal_distribution<double> distribution(-10.0, 10.0);
  for ( int ii=0; ii<NUM; ++ii )
  {
    sus[ii]= RodShapedBacterium(
      Vec3(
        distribution(generator),
        distribution(generator),
        distribution(generator)
      )
    );
    inf[ii]= GenericInfected(
      Vec3(
        distribution(generator),
        distribution(generator),
        distribution(generator)
      )
    );
    phg[ii]= Phage(
      Vec3(
        distribution(generator),
        distribution(generator),
        distribution(generator)
      )
    );
  }
}

void testGrid()
{
  tg.createCLL(sus,inf,phg);
  sus_check=sus;
  // testAllNieghbours();
}

int main(int argc, char const *argv[])
{
  std::cout << BOOST_LIB_VERSION << '\n';
  setUp();
  testGrid();
  return 0;
}
