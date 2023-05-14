/*
  Script to benchmark vector class performance
  g++ -c MathUtility.cpp -O3 -march=native -std=c++17 &&
  g++ MathUtility.o vec_tests.cpp -o vec_test.out -O3 -march=native -std=c++17

*/

#include <iostream>
#include <cmath>
#include <random>
#include <vector>

#include "timer.hpp"

#include "MathUtility.hpp"
// #include "RodShapedBacteria.hpp"
// #include "biofilm_class.hpp"
#include "test.hpp"
#include "constants.hpp"

namespace Vec3Ops
{
  template <class T>
  T addVec3(const T& v_1, const T& v_2)
  {
    return
    T{
      v_1[0] + v_2[0],
      v_1[1] + v_2[1],
      v_1[2] + v_2[2]
    };
  }

  template <class T>
  T multVec3(const T& v_1, const T& v_2)
  {
    return
    T{
      v_1[0] * v_2[0],
      v_1[1] * v_2[1],
      v_1[2] * v_2[2]
    };
  }

  template <class T>
  double sumVec3(const T& v_1)
  {
    return {v_1[0]+v_1[1]+v_1[2]};
  }

  template <class T>
  double dotVec3(const T &v1, const T &v2)
  {
    return sumVec3( v1 * v2 );
  }

  template <class T>
  T crossVec3(const T &v1, const T &v2)
  {
    return
    {
      v1[1] * v2[2] - v1[2] * v2[1],
      v1[2] * v2[0] - v1[0] * v2[2],
      v1[0] * v2[1] - v1[1] * v2[0]
    };
  }
}

int main(int argc, char const *argv[]) {
  /*
    Generate a list of 10000 random vectors for trialing. This will then be
    compared to the eqivalent numbers for random vectors.
  */
  std::mt19937 generator { 32536 };
  std::uniform_real_distribution<double> rndm_dist(-10000,10000);

  constexpr long test_its { 2064 };

  std::vector<Vec3> v1; v1.reserve(test_its);
  std::vector<Vec3> v2; v2.reserve(test_its);
  std::vector<std::array<double, 3> > tst_vec1; tst_vec1.reserve(test_its);
  std::vector<std::array<double, 3> > tst_vec2; tst_vec2.reserve(test_its);

  for (int ii=0; ii<test_its;++ii)
  {
    v1.emplace_back(
      rndm_dist(generator),
      rndm_dist(generator),
      rndm_dist(generator)
    );
    v2.emplace_back(
      rndm_dist(generator),
      rndm_dist(generator),
      rndm_dist(generator)
    );
    for ( int jj=0; jj<3; ++jj)
    {
      tst_vec1[ii][jj]=v1[ii][jj];
      tst_vec2[ii][jj]=v2[ii][jj];
    }
  }

  // std::array<double, 3> tmp_stl_arr;
  // Timer t_stl;
  // for (int ii=0; ii<test_its; ++ii)
  // {
  //   tmp_stl_arr=addVec3(tst_vec1[ii],tst_vec2[ii]);
  // }
  // t_stl.accumulate();
  // std::cout << tmp_stl_arr[0] << " "
  //           << tmp_stl_arr[1] << " "
  //           << tmp_stl_arr[2] << " "
  //           << '\n';
  // std::cout << "stl add: " << t_stl << '\n';

  double stl_dot;
  Timer t_stl_dot;
  for (int ii=0; ii<test_its; ++ii)
  {
    stl_dot+=Vec3Ops::dotVec3(tst_vec1[ii],tst_vec2[ii]);
  }
  t_stl_dot.accumulate();
  std::cout << stl_dot
            << '\n';
  std::cout << "stl dot: " << t_stl_dot << '\n';

  // Vec3 tmp;
  // Timer t_add;
  // for (int ii=0; ii<test_its; ++ii)
  // {
  //   tmp=addVec3(v1[ii],v2[ii]);
  // }
  // t_add.accumulate();
  // std::cout << "tmp: " << tmp << '\n';
  //
  // std::cout << "vec3 add: " << t_add << '\n';

  double vec3;
  Timer t_vec3;
  for (int ii=0; ii<test_its; ++ii)
  {
    vec3+=Vec3Ops::dotVec3(tst_vec1[ii],tst_vec2[ii]);
  }
  t_vec3.accumulate();
  std::cout << vec3
            << '\n';
  std::cout << "vec3 dot: " << t_vec3 << '\n';

  // Biofilm<RodShapedBacterium> filmC{};
  // std::string filename_long { "/media/rory/Elements/biofilm_37554/output/data/vis_biofilm_00455.txt" };
  // RodShapedBacterium::setMTSeed( 234345 );
  // filmC.populateCellsFromFile(filename_long);
  // Timer t_c;
  // // filmC.getAllForceTorque();
  // std::cout << "time: " << t_c.elapsed() << " for "
  //           << filmC.getCellsInFilm().size()
  //           << " cells" << '\n';
  //
  // Timer t_cells;
  // stl_dot=0;
  // std::vector<RodShapedBacterium>& cells { filmC.getCellsInFilm() };
  // std::vector< std::array<double,3> > c; c.reserve(test_its);
  // for ( int ii=0; ii < test_its; ++ii )
  // {
  //   Vec3 ll { cells[ii].getRcm() };
  //   for ( int jj=0; jj<3; ++jj)
  //   {
  //     c[ii][jj]=ll[jj];
  //   }
  // }
  //
  //
  // for ( int ii=0; ii < test_its; ++ii )
  // {
  //   stl_dot+=dot(c[ii],c[ii]);
  // }
  // t_cells.accumulate();
  // std::cout << stl_dot
  //           << '\n';
  // std::cout << "cells dot: " << t_cells << '\n';

  return 0;
}
