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

// TODO: Make parallel versions of the neighbour list functions

#ifndef GRID_CLASS
#define GRID_CLASS

// Standard libraries
#include <iostream>
#include <cmath>
#include <vector>
#include <array>
#include <random>
#include <memory>
#include <algorithm>

// Custom libraries
#include "constants.hpp"
#include "utilities.hpp"
#include "MathUtility.hpp"

template <typename P>
void getMinMax(
  const std::vector<P*> &pars,
  Vec3   &min_pos,
  Vec3   &max_pos,
  double &max_size
)
{
  max_size= -constants::INF;
  min_pos = Vec3 { constants::INF, constants::INF, constants::INF};
  max_pos = -min_pos;
  for ( auto pp : pars )
  {
    min_pos  = min3(min_pos,pp->getPos());
    max_pos  = max3(max_pos,pp->getPos());
    max_size = std::max(max_size,2*pp->getEffectiveR());
  }
}

/**
  @brief
  Coarse grid for use in cell linked list creation
*/
template < class S >
class GridCells
{

public:
  Uint3 mGridSize { 0,0,0 };            //!< Number of elements in grid
  uint  mNumGridCells { 0 };            //!< Total number of grid cells
  Vec3  mCellSize { 0,0,0 };            //!< Side lengths of the grid cell
  Vec3  mOrigin { 0,0,0 };              //!< Grid coordinate origin
  double mVerletRadius { 0 };           //!< Radius to be a neighbour
  double mVerletSkin { 0 };             //!< Skin length extra to the interaction length
  double mVerletMoveR { 0 };            //!< Move more than this, redo NN lists
  uint  mVerletUptime{ 100000000 };     //!< Rebin after mVerletUptime timesteps

  uint getParticleCellID(const Uint3 grid_cell)
  {
    return mGridSize.y*mGridSize.x*grid_cell.z +
           mGridSize.x*grid_cell.y +
           grid_cell.x;
  }

  Uint3 getCellIndex(const Vec3 v)
  {
    return
    Uint3(
      floor( static_cast <uint> ( ( v.x - mOrigin.x ) / mCellSize.x ) ),
      floor( static_cast <uint> ( ( v.y - mOrigin.y ) / mCellSize.y ) ),
      floor( static_cast <uint> ( ( v.z - mOrigin.z ) / mCellSize.z ) )
    );
  }

  struct GridCell {
    S* mHead { nullptr };             //!< First cell in this GridCell
  };
  std::vector< GridCell > mGridCells;//!< GridCells in this grid

  GridCells(std::vector<S*> &cells,double _cell_size)
  {
    // mCellSize.x=mCellSize.y=mCellSize.z=S::mAvgDivLen+2*S::mRadius;
    mCellSize.x   = mCellSize.y = mCellSize.z = _cell_size;
    mVerletSkin   = 0.1*constants::nondim_rodSpheroDiam;
    mVerletRadius = mVerletSkin+_cell_size;
    mVerletMoveR  = 0.5*mVerletSkin;
    allocateGrid(cells);
  }

  /**
    Parameters:
      particle_list: Vector of pointers to all the particles to add to the grid cells
    Effect:
      Create a double linked list to create chains of particles associated with
      each grid node
  */
  void updateGridCellLists(
    const std::vector< S* > &particle_list
  )
  {
    for ( auto &grid_cell : mGridCells )
    {
      grid_cell.mHead = nullptr;
    }

    for ( auto &particle : particle_list )
    {
      uint grid_index { getParticleCellID( getCellIndex(particle->getPos()) ) };

      particle->setHeadLink(mGridCells[grid_index].mHead);

      particle->getGridCell() = grid_index;

      // The grid node always points to the last added particle
      mGridCells[grid_index].mHead = particle;

    }
  }

  void createNeighbourLists(
    const std::vector< S* > &particle_list
  )
  {
    for ( auto &particle : particle_list )
    {
      particle->getNeighbourList().clear(); particle->getNeighbourList().reserve(6);

      // Get the grid node this particle lives in

      Uint3 cell_indices { getCellIndex(particle->getPos()) } ;

      // loop over all the neighbour particles
      for ( int iz=umax(cell_indices.z-1,0); iz<=umin(cell_indices.z+1,mGridSize.z-1); ++iz )
      {
        for ( int iy=umax(cell_indices.y-1,0); iy<=umin(cell_indices.y+1,mGridSize.y-1); ++iy )
        {
          for ( int ix=umax(cell_indices.x-1,0); ix<=umin(cell_indices.x+1,mGridSize.x-1); ++ix )
          {
            Uint3 neighbour_grid_coords = Uint3(ix,iy,iz);

            // Get the first particle in this node
            S* other_particle {
              mGridCells[getParticleCellID(neighbour_grid_coords)].mHead
            };
            while ( other_particle != nullptr ) {
              // Same particle is matching type and id
              bool same_particle {
                ( particle->getID() == other_particle->getID() )
                &&
                ( particle->getMyType() == other_particle->getMyType() )
              };
              const Vec3 sep { particle->getPos()-other_particle->getPos() };

              const double max_sep {
                particle->getEffectiveR()+
                other_particle->getEffectiveR()+
                mVerletSkin
              };
              const bool in_range {
                dot2( sep ) <= dot2( max_sep )
              };

              if ( in_range && !( same_particle ) )
              {
                particle->getNeighbourList().push_back( other_particle );
              }
              // Go down the chain of linked particles in the node
              other_particle = other_particle->getHeadLink();
            }
          }
        }
      }
    }
  }

  void allocateGrid(
    const std::vector< S* > &particle_list
  )
  {
    Vec3 min_pos,max_pos;
    double max_size;
    getMinMax(particle_list,min_pos,max_pos,max_size);

    assert(
      (max_size<=mCellSize.x) && (max_size<=mCellSize.y) && (max_size<=mCellSize.z)
    );

    mGridSize = Uint3(
      findNextPowerOf2(floor( 1 + 2*std::max(fabs(max_pos.x),fabs(min_pos.x)) / mCellSize.x )),
      findNextPowerOf2(floor( 1 + 2*std::max(fabs(max_pos.y),fabs(min_pos.y)) / mCellSize.y )),
      findNextPowerOf2(floor( 1 + 2*std::max(fabs(max_pos.z),fabs(min_pos.z)) / mCellSize.z ))
    );
    mOrigin=Vec3(
      -0.5*mGridSize.x*mCellSize.x,
      -0.5*mGridSize.y*mCellSize.y,
      -0.5*mGridSize.z*mCellSize.z
    );
    for ( auto &gg : mGridCells )
    {
      gg.mHead = nullptr;
    }
    mNumGridCells = mGridSize.x*mGridSize.y*mGridSize.z;
    mGridCells.resize(mNumGridCells);
  }

  void updateVerletList(
    const std::vector< S* > &particle_list
  )
  {
    allocateGrid(particle_list);

    updateGridCellLists(particle_list);

    createNeighbourLists(particle_list);

    // If a new list is created, store the cell's current positions
    for ( auto &particle : particle_list )
      particle->setLoggedPos();
  }

  ~GridCells()
  {}
};

#endif
