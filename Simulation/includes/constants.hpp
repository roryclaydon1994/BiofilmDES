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

#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <cmath>

// Note: usually modE is (E/(1-nu^2)) but assuming incompressible bacteria so
// nu=1/2. This cancels another factor of 4/3 in the Hertzian force calc

// Note: time units use the rod average modE and friction per unit length atm
// The length scale is also set using the rod average diameter (1 micron)

namespace constants
{

  // --------------------------- Common To All -----------------------------------
  inline constexpr double zeta{ 200 }; // drag coefficient Pa h
  inline constexpr double effective_g { 0.01 };
  inline constexpr double baseModE { 4e6 };            // Pa
  inline constexpr double baseLengthScale{ 1 };        // 1 micron
  inline constexpr double surfaceModE{ 100*baseModE }; // Pa (agar effective stiffness)

  inline constexpr double nondim_min_z { -0.5 / baseLengthScale }; // location of agar
  inline constexpr double nondim_agar_roughness { 0.01 / baseLengthScale }; // location of agar
  inline constexpr double nondim_agar_mod_E { surfaceModE / baseModE };
  inline constexpr double baseTimeScale { zeta / baseModE };       // hours
  inline constexpr double dt { 5e-7 }; // hours
  inline constexpr double nondim_dt { dt / baseTimeScale };

  // Mathematical constants
  inline constexpr double pi{3.14159265358979323846};

  // Random number generator seed
  inline constexpr long SEED { 31415926535 };

  // Maximum allowed value for a uint to take in this Simulation
  inline constexpr uint INF { 0xffffffff };

  // ------------------------ RodShapedBacterium -------------------------------

  // Dimensional
  inline constexpr double rodGrwthRtePreFac{ 4.0 };         // microns / hour
  inline constexpr double rodModE{ baseModE };              // Pa
  inline constexpr double rodSpheroDiam{ baseLengthScale }; // microns
  inline constexpr double avg_div_L{ 5.0 };                 // microns

  // Non-dimensional
  inline constexpr double nondim_rodGrwthRtePreFac
  {
    rodGrwthRtePreFac * baseTimeScale / baseLengthScale
  };
  inline constexpr double nondim_rodModE{ rodModE / baseModE };
  inline constexpr double nondim_rodSpheroDiam{
    rodSpheroDiam / baseLengthScale
  };
  inline constexpr double nondim_avg_div_L{ avg_div_L / baseLengthScale };
  inline constexpr double nondim_init_length {
    0.5 * ( nondim_avg_div_L - nondim_rodSpheroDiam )
  };

  // ------------------------ SphericalBacterium -----------------------------
  inline constexpr double sphericalModE{ 10e6 };             // Pa
  inline constexpr double spherical_init_radius{ 0.782/pow(2.0,1.0/3.0) }; // microns
  inline constexpr double avg_div_radius{ 0.782 }; // microns
  inline constexpr double sphericalGrwthRtePreFac{
    2*(avg_div_radius-spherical_init_radius)
  }; // microns / hour

  // RodShapedBacterium cell specific constants non-dimensionalised
  inline constexpr double nondim_sphericalGrwthRtePreFac
  {
    sphericalGrwthRtePreFac * baseTimeScale / baseLengthScale
  };
  inline constexpr double nondim_sphericalModE{ sphericalModE/baseModE };
  inline constexpr double nondim_spherical_init_radius{
    spherical_init_radius/baseLengthScale
  };
  inline constexpr double nondim_avg_div_radius{
    avg_div_radius / baseLengthScale
  };

  // ------------------------------- Phage -------------------------------------
// #ifdef PHAGE
  inline constexpr double phage_diffusion { 4.0 }; // m^2/s

  inline constexpr double nondim_phage_diffusion {
    ( phage_diffusion*3600.0 ) * baseTimeScale / pow(baseLengthScale,2)
  };

  // https://www.ncbi.nlm.nih.gov/pmc/articles/PMC93561/pdf/jb001677.pdf
  inline constexpr double lysis_period { 23.6 }; // minutes
  inline constexpr double std_lysis_period { 4.3 }; // minutes
  inline constexpr double avg_eclipse_period { 18.3 }; // minutes
  inline constexpr double std_eclipse_period {  3.4 }; // minutes

  inline constexpr double nondim_lysis_period {
    ( lysis_period/60.0 ) / baseTimeScale
  };
  inline constexpr double nondim_std_lysis_period {
    ( std_lysis_period/60.0 ) / baseTimeScale
  };

  inline constexpr double nondim_avg_eclipse_period {
    ( avg_eclipse_period/60.0 ) / baseTimeScale
  };
  inline constexpr double nondim_std_eclipse_period {
    ( std_eclipse_period/60.0 ) / baseTimeScale
  };

  inline constexpr double nondim_avg_lysis_rate { 1.0 / nondim_lysis_period };

  inline constexpr uint burst_size { 50 }; // not used if variable burst used
  inline constexpr double phage_production_rate =
    247.1/( nondim_lysis_period-nondim_avg_eclipse_period );

// #endif // end phage definition
  // ------------------------------ Adhesion -----------------------------------

#ifdef ADHESION
  // Dimensional
  inline constexpr double EPS_radius { 0.5 };                        // microns
  inline constexpr double Rd { 0.5*rodSpheroDiam+EPS_radius };       // microns
  inline constexpr double Ri { Rd*2 };                               // microns

  // Non dimensional
  inline constexpr double nondim_Rd { Rd / baseLengthScale };
  inline constexpr double nondim_Ri { Ri / baseLengthScale };
  inline constexpr double kappa_depletion { 5e-3 };
  inline constexpr double repulsion_strength { 10 };

#endif // end adhesion definition
  // -------------------------- Chaining & AG43 -------------------------------

  // Assuming a perfectly elastic rod linking the two RodShapedBacterium, with a Young's
  // modulus E_rod = 5300 Pa, radius = 0.5 microns, with a circular cross-section
  // the bending stiffness is given by K = R_rod (pi r^2 / 2)^2
  // inline constexpr double nondim_K_bend {
  //   5300*pi*0.25*(0.1*0.1*0.1*0.1) / (rodModE)
  // };
  // inline constexpr double kappa { 1 };
  inline constexpr double nondim_K_bend {
    1
  };
  // inline constexpr double nondim_kappa { kappa / ( baseModE * rodSpheroDiam ) };
  inline constexpr double nondim_kappa { 1 };

  inline constexpr double max_extension { 0.1*nondim_rodSpheroDiam };

  // ------------------------- Integration Constants ---------------------------

  // Biofilm specific constants
  inline constexpr double colony_size_target{ 37500.0 };
  // inline constexpr double colony_size_target{ 30.0 };
  // inline constexpr double colony_size_target{ 1200.0 };

  // Grid
  inline constexpr double box_width { nondim_rodSpheroDiam+nondim_avg_div_L };
  // periodic boundary for phage
  inline constexpr double max_grid_size { 200*box_width };

}
#endif // End fileguard
