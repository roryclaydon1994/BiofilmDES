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
#include <array>
#include <cassert>
#include <initializer_list>
#include <iostream>

struct Uint3
{
  uint x,y,z;
  Uint3(uint _x, uint _y, uint _z) : x{_x},y{_y},z{_z} {}

  Uint3& operator+= (const Uint3 &v1);
  Uint3& operator-= (const Uint3 &v1);
  Uint3& operator*= (const uint v1);

  // Display vector
  friend std::ostream& operator<< (std::ostream &out, const Uint3 &vec);
};
Uint3 operator-(const Uint3 &v1);
Uint3 operator+(const Uint3 &v1, const Uint3 &v2);
Uint3 operator-(const Uint3 &v1, const Uint3 &v2);
Uint3 operator*(const Uint3 &v1, const uint scal);
Uint3 operator*(const uint scal, const Uint3 &v1);

struct Int3
{
  int x,y,z;
  Int3(int _x, int _y, int _z) : x{_x},y{_y},z{_z} {}

  Int3& operator+= (const Int3 &v1);
  Int3& operator-= (const Int3 &v1);
  Int3& operator*= (const int v1);

  // Display vector
  friend std::ostream& operator<< (std::ostream &out, const Int3 &vec);
};
Int3 operator-(const Int3 &v1);
Int3 operator+(const Int3 &v1, const Int3 &v2);
Int3 operator+(const Int3 &v1, const Uint3 &v2);
Int3 operator+(const Uint3 &v1, const Int3 &v2);
Int3 operator-(const Int3 &v1, const Int3 &v2);
Int3 operator*(const Int3 &v1, const int scal);
Int3 operator*(const int scal, const Int3 &v1);

struct Vec2
{
  double x,y;
  Vec2(double _x, double _y) : x{_x},y{_y} {}

  Vec2& operator+= (const Vec2 &v1);
  Vec2& operator-= (const Vec2 &v1);
  Vec2& operator*= (const double v1);

  // Display vector
  friend std::ostream& operator<< (std::ostream &out, const Vec2 &vec);
};
Vec2 operator-(const Vec2 &v1);
Vec2 operator+(const Vec2 &v1, const Vec2 &v2);
Vec2 operator-(const Vec2 &v1, const Vec2 &v2);
Vec2 operator*(const Vec2 &v1, const double scal);
Vec2 operator*(const double scal, const Vec2 &v1);

class Vec3 {
public:
  // std::array<double,3> mVec;
  double x, y, z;

  Vec3 () : x{0},y{0},z{0} {}
  Vec3 (double x, double y, double z) : x{x},y{y},z{z} {}
  Vec3 (double x) : Vec3{x,x,x} {}
  // Vec3 (std::array<double,3> init_v);   // Copy in an initialiation vector

  // Vec3 (std::initializer_list<double> list); // Brace initialiation

  // Subscripting operator
  // double& operator[](int ii);
  // const double& operator[](int ii) const;

  // += and -= operator
  Vec3& operator+= (const Vec3 &v1);
  Vec3& operator-= (const Vec3 &v1);

  // *= operator
  Vec3& operator*= (const double v1);

  // Get vector
  const std::array<double,3>& getVec() const;

  // Display vector
  friend std::ostream& operator<< (std::ostream &out, const Vec3 &vec);

  // Member vector operators
  double norm() const;

  // Return number of elements of Vec3
  int size() const;

  void zero() { x=y=z=0; }

  virtual ~Vec3 ();
};

/* Comparison operators */
bool operator== (const Vec3 &v1, const Vec3 &v2);
bool operator!= (const Vec3 &v1, const Vec3 &v2);

/* Overloaded arithmetic operators */
Vec3 operator-(const Vec3 &v1);
Vec3 operator+(const Vec3 &v1, const Vec3 &v2);
Vec3 operator-(const Vec3 &v1, const Vec3 &v2);
Vec3 operator*(const Vec3 &v1, const double scal);
Vec3 operator*(const double scal, const Vec3 &v1);
Vec3 operator/(const Vec3 &v1, const double scal);

/* Vector operators between any two instances of Vec3 */
Vec3 min3(Vec3 a, Vec3 b);
Vec3 max3(Vec3 a, Vec3 b);
Vec3 mod3(Vec3 r, double a);

// Find the smallest angle between two lines with direction vectors v1 and v2
double norm(const std::array<double,3> &v);
double norm(const Vec3 &v);
double angularSep(const Vec3 &v1, const Vec3 &v2);
double getThetaFromOrientation(const Vec3 &orientation);

/* This is the definition of alpha in Pournin referenced in biofilm.hpp */
template <typename T> int sgn(T val); // Return the sign of any type T, or 0
double clampRoot(double root);

double getPointToLineSegDist(const Vec3 &p, const Vec3 &r1, const Vec3 &r2);
/**<
  Find the shortest distance from a point to a line segment in d dims

  Based on Dan Sunday's implementation
  http://geomalgorithms.com/a02-_lines.html

  @param[in] p: the location vector of the point
  @param[in] r1: the start of the line segment
  @param[in] r2: the end of the line segment
  @returns dist: shortest distance to the line segment

*/

void closestApproachLineSegmentsPournin(const Vec3 &x_1,
                                        const Vec3 &a_1,
                                        const Vec3 &x_2,
                                        const Vec3 &a_2,
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

template < typename T >
inline double dot(T v1, T v2)
{
  return
  {
    v1.x * v2.x +
    v1.y * v2.y +
    v1.z * v2.z
  };
}

template < typename T >
inline double dot2(T v)
{
  return dot(v,v);
}

template <  >
inline double dot2(double v)
{
  return v*v;
}

template < typename T >
inline T normalize(T v)
{
  return v / sqrt(dot2(v));
}

template < typename T >
inline T cross(T v1, T v2)
{
  return
  {
    v1.y * v2.z - v1.z * v2.y,
    v1.z * v2.x - v1.x * v2.z,
    v1.x * v2.y - v1.y * v2.x
  };
}

inline Vec3 getOrientationFromAngles(Vec3 angles)
{
  return
  {
    sin(angles.y)*cos(angles.x),
    sin(angles.y)*sin(angles.x),
    cos(angles.y)
  };
}

inline Vec2 getAnglesFromOrientation(Vec3 orientation)
{
  return { acos(orientation.z), atan2(orientation.y,orientation.x)  };
}

inline Vec3 projectAngVelRod(const Vec3 ang_vel, const Vec3 angles)
{
  return {
    ang_vel.z,                                          // dphi / dt
    -sin(angles.x)*ang_vel.x + cos(angles.x)*ang_vel.y, // dtheta /dt
    dot(ang_vel,getOrientationFromAngles(angles))       // dpsi / dt
  };
}

inline int umax(int a, int b)
{
  return (a<b) ? b : a;
}

inline int umin(int a, int b)
{
  return (a<b) ? a : b;
}

inline bool isclose(
  const double x,
  const double y,
  const double rel_tol=1e-5,
  const double abs_tol=1e-8
)
{
  // python is close
  return fabs(x-y)<=std::max( rel_tol*std::max(fabs(x),fabs(y)),abs_tol );
}

inline bool isclose(
  const Vec3 a,
  const Vec3 b,
  const double rel_tol=1e-5,
  const double abs_tol=1e-8
)
{
  // python is close
  return isclose(a.x,b.x,rel_tol,abs_tol) &&
         isclose(a.y,b.y,rel_tol,abs_tol) &&
         isclose(a.z,b.z,rel_tol,abs_tol);
}

#endif // end fileguard
