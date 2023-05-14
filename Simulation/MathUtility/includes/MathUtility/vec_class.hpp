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

#ifndef VEC_CLASS
#define VEC_CLASS
/*
  Useful operations on vectors

  Pournin, L., Weber, M., Tsukahara, M. et al.
  Three-dimensional distinct element simulation of spherocylinder crystallization.
  Granular Matter 7, 119–126 (2005). https://doi.org/10.1007/s10035-004-0188-4
  Lumelsky, V.J.:
  On fast computation of distance between line segments.
  Information Processing Letters21, 1985
*/

// Standard libraries
#include <vector>
#include <cassert>
#include <initializer_list>
#include <iostream>

class Rvec {
private:
  std::vector<double> mVec;

public:

  Rvec() = default; // Create basic constructor

  Rvec (int n_elems); // Zero initialiation of vector with n_elements

  Rvec (std::vector<double> init_v);   // Copy in an initialiation vector

  Rvec (std::initializer_list<double> list); // Brace initialiation

  // Subscripting operator
  double& operator[](int ii);
  const double& operator[](int ii) const;

  // += and -= operator
  Rvec& operator+= (const Rvec &v1);
  Rvec& operator-= (const Rvec &v1);

  // *= operator
  Rvec& operator*= (const double v1);

  // Get vector
  std::vector<double> getVec() const;

  // Display vector
  friend std::ostream& operator<< (std::ostream &out, const Rvec &vec);

  // Member vector operators
  double dot(const Rvec &other_vec) const;
  Rvec cross(const Rvec &other_vec) const;
  double norm() const;

  // Return number of elements of Rvec
  int size() const;

  virtual ~Rvec ();
};

/* Comparison operators */
bool operator== (const Rvec &v1, const Rvec &v2);
bool operator!= (const Rvec &v1, const Rvec &v2);

/* Overloaded arithmetic operators */
Rvec operator- (Rvec v1);
Rvec operator+(const Rvec &v1, const Rvec &v2);
Rvec operator-(const Rvec &v1, const Rvec &v2);
Rvec operator*(const Rvec &v1, const double scal);
Rvec operator*(const double scal, const Rvec &v1);
Rvec operator/(const Rvec &v1, const double scal);

/* Vector operators between any two instances of Rvec */
double dot(const Rvec &v1, const Rvec &v2);
Rvec cross(const Rvec &v1, const Rvec &v2);

// Find the smallest angle between two lines with direction vectors v1 and v2
double norm(const std::vector<double> &v);
double norm(const Rvec &v);
double angularSep(const Rvec &v1, const Rvec &v2);
double getThetaFromOrientation(const Rvec &orientation);

/* This is the definition of alpha in Pournin referenced in biofilm.hpp */
template <typename T> int sgn(T val); // Return the sign of any type T, or 0
double clampRoot(double root);
void closestApproachLineSegmentsPournin(const Rvec &x_1,
                                        const Rvec &a_1,
                                        const Rvec &x_2,
                                        const Rvec &a_2,
                                        double &s_star, double &t_star);
/*
  Find the closest distance between two line segments space. This algorithm is
  based off of https://core.ac.uk/download/pdf/147910415.pdf
  Line segment 1 is: r_1 = x_1 + s a_1
  Line segment 2 is: r_2 = x_2 + t a_2
  The length of segment i is L_i, and it is assumed that |a_i| = L_i/2
  Parameters:
    x_1,x_2: The mid point of the line segments
    a_1,a_2: The direction vector of the segment with |a_i| = L_i/2
    s_star, t_star: references to the parameterisation of the two lines resp.

  If the lines are parallel then the proximate points are set to be the
  midpoint of the overlap.

  From [2], thresholding the denominator should be enough to obtain accurate
  results.
*/

#endif // end fileguard
