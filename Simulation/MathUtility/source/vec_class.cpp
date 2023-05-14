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

// Standard libraries
#include <iostream>
#include <cmath>
#include <cassert>

// User defined
#include "vec_class.hpp"

Rvec::Rvec (int n_elems)
{
  assert(n_elems>0);
  mVec.resize(n_elems,0.0);
}

Rvec::Rvec (std::vector<double> init_v)
{
  assert(init_v.empty()==false);
  mVec = init_v;
}

Rvec::Rvec (std::initializer_list<double> list) :
  Rvec(static_cast<int>(list.size()))
{
  int count{ 0 };
  for (auto element : list)
  {
    mVec[count] = element;
    ++count;
  }
}

int Rvec::size() const
{
  return mVec.size();
}

double& Rvec::operator[](int ii)
{
  assert( ii>=0 && ii<this->size());
  return mVec[ii];
}

const double& Rvec::operator[](int ii) const
{
  assert( ii>=0 && ii<this->size());
  return mVec[ii];
}

Rvec& Rvec::operator+= (const Rvec &v1)
{
  Rvec result{ *this + v1};
  *this = result;
  return *this;
}

Rvec& Rvec::operator-= (const Rvec &v1)
{
  return *this+=(-v1);
}

Rvec& Rvec::operator*= (const double v1)
{
  Rvec result{ *this * v1};
  *this = result;
  return *this;
}

std::vector<double> Rvec::getVec() const
{
  return mVec;
}

Rvec operator+(const Rvec &v1, const Rvec &v2)
{
  assert(v1.size() == v2.size());
  Rvec result (static_cast<int>(v1.size()));
  for(int ii = 0; ii < v1.size(); ++ii)
  {
    result[ii] = v1[ii] + v2[ii];
  }
  return result;
}

Rvec operator-(const Rvec &v1, const Rvec &v2)
{
  return Rvec{v1 + (-v2)};
}

Rvec operator*(const Rvec &v1, const double scal)
{
  Rvec res( v1.size() );
  for(int ii = 0; ii < v1.size(); ++ii)
  {
    res[ii] = v1[ii] * scal;
  }
  return res;
}

Rvec operator*(const double scal, const Rvec &v1)
{
  return v1*scal;
}

Rvec operator/(const Rvec &v1, const double scal)
{
  assert(scal != 0); // Avpid division by 0
  return v1 * (1/scal);
}

Rvec::~Rvec()
{}

Rvec operator- (Rvec v1)
{
  Rvec minus_vec(static_cast<int>(v1.size()));
  for (int ii=0; ii < minus_vec.size(); ++ii)
  {
    minus_vec[ii] = -v1[ii];
  }
  return minus_vec;
}

double dot(const Rvec &v1, const Rvec &v2)
{
  assert(v1.size() == v2.size());
  double dot_product{0.0};
  for (int ii=0; ii<v1.size();++ii)
  {
    dot_product += v1[ii] * v2[ii];
  }
  return dot_product;
}

double Rvec::dot(const Rvec &other_vec) const
{
  return dot(*this,other_vec);
}

Rvec cross(const Rvec &v1, const Rvec &v2)
{
  int len_v{v1.size()};
  assert(len_v == v2.size());
  assert(len_v == 3);
  Rvec cross_product(3);

  cross_product[0] = v1[1] * v2[2] - v1[2] * v2[1];
  cross_product[1] = v1[2] * v2[0] - v1[0] * v2[2];
  cross_product[2] = v1[0] * v2[1] - v1[1] * v2[0];

  return cross_product;
}

Rvec Rvec::cross(const Rvec &other_vec) const
{
  return cross(*this,other_vec);
}

std::ostream& operator<< (std::ostream &out, const Rvec &vec)
{
  out << "( ";
  for( auto vm : vec.mVec )
  {
    out << vm << ", ";
  }
  out << ")";
  return out;
}

bool operator== (const Rvec &v1, const Rvec &v2)
{
  if (v1.size()!=v2.size())
  {
    return false;
  }
  else
  {
    for( int ii=0; ii<v1.size(); ++ii )
    {
      if (v1[ii] != v2[ii]) {
        return false;
      }
    }
    return true;
  }
}

bool operator!= (const Rvec &v1, const Rvec &v2)
{
  return !(v1 == v2);
}

double norm(const std::vector<double> &v)
{
  double norm{0.0};
  for ( auto elem : v )
  {
    norm += elem*elem;
  }
  return sqrt(norm);
}

double norm(const Rvec &v)
{
  return norm(v.getVec());
}

double Rvec::norm() const
{
  return norm(mVec);
}

double angularSep(const Rvec &v1, const Rvec &v2)
{
  double scalar_product = v1.dot(v2) / (v1.norm() * v2.norm());
  double theta{acos(scalar_product)};
  return theta;
}

double getThetaFromOrientation(const Rvec &orientation)
{
  return atan2(orientation[1],orientation[0]);
}


// template signum function
template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

// Return the sign of a double x - This is intended ony for use in this file
inline int sgnDouble(double x)
{
  // Return a value with the magnitude 1 and sign of x
  // return std::copysign(1.0,x);
  return sgn(x);
}

double clampRoot(double root)
{
  if (root>=1) return 1;
  else if (root<=-1) return -1;
  else return root;
}

void closestApproachLineSegmentsPournin(const Rvec &x_1,
                                        const Rvec &a_1,
                                        const Rvec &x_2,
                                        const Rvec &a_2,
                                        double &s_star, double &t_star)
{
  /*
    Minimise the distance between two line segments.
    Define R(s,t) = a s^2 + 2 b s t + c t^2 + 2 d s + 2 e t + f
                  = |x_1 - x_2|^2
    Return the values s^* and t^* s.t. R(s^*,t^*) is minimal subject to the
    constraint that s,t \in [-1,1]^2

    If the lines are parallel then the proximate points are set to be the
    midpoint of the overlap.

    Algorithm from [1] and [2]
  */
  Rvec x_21 = x_1 - x_2; // vector from x2 to x1
  double a =  dot(a_1,a_1);
  double b = -dot(a_1,a_2);
  double c =  dot(a_2,a_2);
  double d =  dot(a_1,x_21);
  double e = -dot(a_2,x_21);
  double delta = a*c - b*b;  // \equiv modulus( a_1 \cross a_2 )

  // Do not consider degenerate cases (one or more segments are points)
  assert(a>0); assert(c>0);
  if ( delta == 0 )  // Parallel
  {
    double mod_a{ sqrt(a) };
    double mod_c{ sqrt(c) };
    double Delta{ (fabs(d) / mod_a) };

    // Define segment 1 end points and centre at 0
    double q_0{-mod_a}; double q_1{mod_a};

    // Define segment 2 end points and centre at sep
    double p_0{Delta - mod_c};
    double p_1{Delta + mod_c};

    if ( q_1 <= p_0 ) {
      // no overlap
      s_star = sgn( -d );
      t_star = sgn( -e );
    }
    else if (q_0 <= p_0 && q_1 >= p_1)
    {
      // segment p is inside q
      // std::cout << "p in q" << '\n';
      t_star = 0;
      s_star = sgn( -d ) * Delta/mod_a;
    }
    else if (q_0 > p_0 && q_1 < p_1)
    {
      // segment q is inside p
      // std::cout << "q in p" << '\n';
      t_star = sgn( -e ) * Delta/mod_c;
      s_star = 0;
    }
    else if (q_0 <= p_0 && q_1 < p_1)
    {
      // skew overlap
      // std::cout << "skew overlap" << '\n';
      double midpoint_of_overlap{ 0.5*( Delta - mod_c + mod_a ) };
      s_star= sgn( -d ) * midpoint_of_overlap/mod_a;
      t_star= sgn( -e ) * (Delta-midpoint_of_overlap)/mod_c;
    }
    else
    {
      std::cout << "FATAL_ERROR! Parallel segments undefined" << '\n';
      exit(2);
    }

    std::cout << "===========================================" << '\n';
    std::cout << "Compare previous parallel method to current" << '\n';
    std::cout << "previous" << '\n';
    std::cout << "s: " << s_star << " t: " << t_star << '\n';
    std::cout << "new" << '\n';
    s_star = 0.5*( clampRoot( (c-e)/b )+clampRoot( -(c+e)/b ) );
    t_star = clampRoot(-(e+b*s_star)/c);
    std::cout << "s: " << s_star << " t: " << t_star << '\n';
    exit(1);
  }
  else // Non-parallel
  {

    t_star = clampRoot( ( b*d - a*e ) / delta );

    double s_star_tmp = (-b*t_star - d) / a;


    s_star = clampRoot( s_star_tmp );

    if ( s_star_tmp < -1 || s_star_tmp > 1 )
      t_star = clampRoot( ( -b*s_star - e ) / c );
  }

}
