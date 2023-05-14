/*
  Templated code for non polymorphic behaviour - depricated.
*/

// -----------------------------------------------------------------------------
// Method using an index grid
// TODO: Need to make the multi indexgrid take arbitrary number of different
// particle vectors
// -----------------------------------------------------------------------------

#ifndef INDEX_GRID_H
#define INDEX_GRID_H

// Standard libraries
#include <iostream>
#include <cmath>
#include <vector>
#include <array>
#include <random>
#include <memory>
#include <algorithm>
//#include <execution>
#include <numeric>
#include <utility>
#include <functional>

// Third party
// #include <boost/compute/algorithm/transform_reduce.hpp>
// #include <boost/compute/algorithm/sort_by_key.hpp>

// User defined
// #include "RodShapedBacteria.hpp"
// #include "chaining_susceptible_class.hpp"
// #include "infected_class.hpp"
// #include "Phage.hpp"
#include "utilities.hpp"
#include "HyperParams.hpp"

extern HyperParams hyperParams;

inline Uint3 getCellIndex(const Vec3 v)
{
  return
  Uint3(
    floor( static_cast <uint> ( ( v.x - hyperParams.origin.x ) / hyperParams.cellSize.x ) ),
    floor( static_cast <uint> ( ( v.y - hyperParams.origin.y ) / hyperParams.cellSize.y ) ),
    floor( static_cast <uint> ( ( v.z - hyperParams.origin.z ) / hyperParams.cellSize.z ) )
  );
}

inline uint getParticleCellID(const Uint3 grid_cell)
{
  // std::cout << "grid index: " << grid_cell << " @ " << v << '\n';
  return hyperParams.gridSize.y*hyperParams.gridSize.x*grid_cell.z +
         hyperParams.gridSize.x*grid_cell.y +
         grid_cell.x;
}

// template < class P >
// void binParticlesToGrid(
//   std::vector<P> &pars
// )
// {
//   std::for_each(
//     std::execution::par_unseq,
//     pars.begin(),
//     pars.end(),
//     [](auto &A)
//     {
//       A.mMyGridCell = getParticleCellID( getCellIndex(A.mPos) );
//     }
//   );
//   // TODO: for each? - Could be removed and put into the lambda
//   // for ( uint ii=0; ii<pars.size(); ++ii )
//   // {
//   //   Uint3 grid_cell { getCellIndex(pars[ii].mPos) };
//   //   pars[ii].mMyGridCell = getParticleCellID(grid_cell);
//   //     // std::cout << "pars["<<ii<<"].mMyGridCell = "
//   //     //           << getParticleCellID(pars[ii].mPos)
//   //     //           << " @ " << pars[ii].mPos << '\n';
//   // }
//   // exit(1);
//   // boost::compute::sort_by_key(cell_list.begin(),cell_list.end(),pars.begin());
//   return;
// }

// template < class P >
// void reorderParticles(
//   std::vector<P> &pars
// )
// {
//   std::sort(
//     std::execution::par_unseq,
//     pars.begin(),
//     pars.end(),
//     [](auto &A, auto &B) { return A.mMyGridCell < B.mMyGridCell; }
//   );
// }
//
// template <>
// inline void reorderParticles<ChainingRodShapedBacterium>(
//   std::vector<ChainingRodShapedBacterium> &pars
// )
// {
//   // std::sort(
//   //   std::execution::par_unseq,
//   //   pars.begin(),
//   //   pars.end(),
//   //   [](auto &A, auto &B) { return A.mMyGridCell < B.mMyGridCell; }
//   // );
// }

template <typename P>
void getMinMax3(std::vector<P> &pars, Vec3 &min_pos, Vec3 &max_pos)
{
  // TODO Replace with Boost/par ex
  for ( auto pp : pars )
  {
    min_pos = min3(min_pos,pp.mPos);
    max_pos = max3(max_pos,pp.mPos);
  }
  // auto tt = std::transform_reduce(
  //   sus.begin(),
  //   sus.end(),
  //   Vec3 { 0,0,0 },
  //   [](S &s) { return s.mPos; }                    // Transform
  //   [](Vec3 a, Vec3 b) { return max3(a,b); }       // Reduce
  // );
}

template < typename P >
void fillInCellBeginEnd(
  std::vector<P> &par,
  std::vector< uint > &cell_begin,
  std::vector< uint > &cell_end
)
{
  uint my_cell, neighbour_cell;
  // std::cout << "pp: " << par.size() << " ng: " << hyperParams.numGridCells << '\n';
  for ( uint id=0; id<par.size(); ++id )
  {
    my_cell = par[id].mMyGridCell;
    // std::cout << "id: " << id << " my_cell: " << my_cell << '\n';

    // This must be the first particle entry
    if ( id==0 )
    {
      // Any choice is safe as if id==0, we never use neighbour_cell
      neighbour_cell = my_cell;
    }
    else
    {
      // Compare to next particle
      neighbour_cell = par[id-1].mMyGridCell;
      // std::cout << "neighbour_cell: " << neighbour_cell << '\n';
    }

    if ( my_cell!=neighbour_cell || id==0 )
    {
      // The previous particle is in a different cell, so this particle is first
      // std::cout << "cell_begin["<< my_cell <<"]=" << id << '\n';
      cell_begin[my_cell]=id;
      // But that means that the previous particle must have been the last
      if ( id > 0 )
      {
        // std::cout << "cell_end["<<neighbour_cell<<"]="<<id << '\n';
        cell_end[neighbour_cell]=id;
      }
    }

    if ( id==par.size()-1 )
    {
      // std::cout << "cell_end["<< my_cell <<"]=" << id+1 << '\n';
      cell_end[my_cell]=id+1;
    }
  }
  return;
}

class IndexGrid
{
public:
  enum Types { SUS, INF, PHG };
  std::array< std::vector< uint >, 3 > mCellBegin;
  std::array< std::vector< uint >, 3 > mCellEnd;

  template < typename S >
  void allocateGrid(
    std::vector<S> &sus
  )
  {
    // TODO: populate from here
    Vec3 min_pos {  constants::INF, constants::INF, constants::INF };
    Vec3 max_pos { -min_pos };
    getMinMax3(sus,min_pos,max_pos);

    // TODO: create a halo and only reallocate grid if a particle enters the halo
    // if ( reallocateGrid(min_pos,max_pos) )
    // std::cout << "cell size: " << hyperParams.cellSize << '\n';
    hyperParams.gridSize = Uint3(
      findNextPowerOf2(floor( 1 + 2*std::max(fabs(max_pos.x),fabs(min_pos.x)) / hyperParams.cellSize.x )),
      findNextPowerOf2(floor( 1 + 2*std::max(fabs(max_pos.y),fabs(min_pos.y)) / hyperParams.cellSize.y )),
      findNextPowerOf2(floor( 1 + 2*std::max(fabs(max_pos.z),fabs(min_pos.z)) / hyperParams.cellSize.z ))
    );
    // std::cout << "hyperParams.gridSize: " << hyperParams.gridSize << '\n';

    hyperParams.numGridCells = hyperParams.gridSize.x*hyperParams.gridSize.y*hyperParams.gridSize.z;
    mCellBegin[SUS].resize(hyperParams.numGridCells);
    std::fill(mCellBegin[SUS].begin(),mCellBegin[SUS].end(),constants::INF);
    mCellEnd[SUS].resize(hyperParams.numGridCells);
    hyperParams.origin=Vec3(
      -0.5*hyperParams.gridSize.x*hyperParams.cellSize.x,
      -0.5*hyperParams.gridSize.y*hyperParams.cellSize.y,
      -0.5*hyperParams.gridSize.z*hyperParams.cellSize.z
    );
    // std::cout << "hyperParams.origin: " << hyperParams.origin << '\n';
    return;
  }

  template < typename S, typename I, typename P >
  void allocateGrid(
    std::vector<S> &sus,
    std::vector<I> &inf,
    std::vector<P> &phg
  )
  {
    // TODO: populate from here
    Vec3 min_pos {  constants::INF, constants::INF, constants::INF };
    Vec3 max_pos { -min_pos };
    getMinMax3(sus,min_pos,max_pos);
    getMinMax3(inf,min_pos,max_pos);
    getMinMax3(phg,min_pos,max_pos);

    // TODO: create a halo and only reallocate grid if a particle enters the halo
    // if ( reallocateGrid(min_pos,max_pos) )
    // std::cout << "cell size: " << hyperParams.cellSize << '\n';
    hyperParams.gridSize = Uint3(
      findNextPowerOf2(floor( 1 + 2*std::max(fabs(max_pos.x),fabs(min_pos.x)) / hyperParams.cellSize.x )),
      findNextPowerOf2(floor( 1 + 2*std::max(fabs(max_pos.y),fabs(min_pos.y)) / hyperParams.cellSize.y )),
      findNextPowerOf2(floor( 1 + 2*std::max(fabs(max_pos.z),fabs(min_pos.z)) / hyperParams.cellSize.z ))
    );
    // std::cout << "hyperParams.gridSize: " << hyperParams.gridSize << '\n';

    hyperParams.numGridCells = hyperParams.gridSize.x*hyperParams.gridSize.y*hyperParams.gridSize.z;
    mCellBegin[SUS].resize(hyperParams.numGridCells);
    mCellBegin[INF].resize(hyperParams.numGridCells);
    mCellBegin[PHG].resize(hyperParams.numGridCells);
    std::fill(mCellBegin[SUS].begin(),mCellBegin[SUS].end(),constants::INF);
    std::fill(mCellBegin[INF].begin(),mCellBegin[INF].end(),constants::INF);
    std::fill(mCellBegin[PHG].begin(),mCellBegin[PHG].end(),constants::INF);
    mCellEnd[SUS].resize(hyperParams.numGridCells);
    mCellEnd[INF].resize(hyperParams.numGridCells);
    mCellEnd[PHG].resize(hyperParams.numGridCells);

    hyperParams.origin=Vec3(
      -0.5*hyperParams.gridSize.x*hyperParams.cellSize.x,
      -0.5*hyperParams.gridSize.y*hyperParams.cellSize.y,
      -0.5*hyperParams.gridSize.z*hyperParams.cellSize.z
    );
    // std::cout << "hyperParams.origin: " << hyperParams.origin << '\n';
    return;
  }


  template < typename S >
  void createCLL(
    std::vector<S> &sus
  )
  {
    allocateGrid(sus);

    binParticlesToGrid(sus);
    reorderParticles(sus);

    fillInCellBeginEnd(sus,mCellBegin[SUS],mCellEnd[SUS]);
    return;
  }

  template < typename S, typename I, typename P >
  void createCLL(
    std::vector<S> &sus,
    std::vector<I> &inf,
    std::vector<P> &phg
  )
  {
    allocateGrid(sus,inf,phg);

    binParticlesToGrid(sus);
    reorderParticles(sus);

    binParticlesToGrid(inf);
    reorderParticles(inf);

    binParticlesToGrid(phg);
    reorderParticles(phg);

    fillInCellBeginEnd(sus,mCellBegin[SUS],mCellEnd[SUS]);
    fillInCellBeginEnd(inf,mCellBegin[INF],mCellEnd[INF]);
    fillInCellBeginEnd(phg,mCellBegin[PHG],mCellEnd[PHG]);
    return;
  }
};

// template void IndexGrid::allocateGrid(
//   std::vector<RodShapedBacterium>     &sus
// );
// template void IndexGrid::allocateGrid(
//   std::vector<RodShapedBacterium>     &sus,
//   std::vector<GenericInfected> &inf,
//   std::vector<Phage>    &phg
// );

// template void IndexGrid::createCLL(
//   std::vector<RodShapedBacterium>     &sus
// );
// template void IndexGrid::createCLL(
//   std::vector<RodShapedBacterium>     &sus,
//   std::vector<GenericInfected> &inf,
//   std::vector<Phage>    &phg
// );

// template void binParticlesToGrid(
//   std::vector<RodShapedBacterium> &pars
// );
// template void reorderParticles(
//   std::vector<RodShapedBacterium> &pars
// );
// template void binParticlesToGrid(
//   std::vector<GenericInfected> &pars
// );
// template void reorderParticles(
//   std::vector<GenericInfected> &pars
// );
// template void binParticlesToGrid(
//   std::vector<Phage> &pars
// );
// template void reorderParticles(
//   std::vector<Phage> &pars
// );
// template void getMinMax3(
//   std::vector<RodShapedBacterium> &pars,
//   Vec3 &min_pos,
//   Vec3 &max_pos
// );
// template void fillInCellBeginEnd(
//   std::vector<RodShapedBacterium> &par,
//   std::vector< uint > &cell_begin,
//   std::vector< uint > &cell_end
// );
// template void fillInCellBeginEnd(
//   std::vector<GenericInfected> &par,
//   std::vector< uint > &cell_begin,
//   std::vector< uint > &cell_end
// );
// template void fillInCellBeginEnd(
//   std::vector<Phage> &par,
//   std::vector< uint > &cell_begin,
//   std::vector< uint > &cell_end
// );

#endif
