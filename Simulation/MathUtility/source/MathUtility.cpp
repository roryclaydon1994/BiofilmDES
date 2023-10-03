// Standard libraries
#include <iostream>
#include <cmath>
#include <cassert>
#include <array>

// User defined
#include "MathUtility.hpp"

Vec2& Vec2::operator+= (const Vec2 &v1)
{
  Vec2 result{ *this + v1};
  *this = result;
  return *this;
}

Vec2& Vec2::operator-= (const Vec2 &v1)
{
  return *this+=(-v1);
}

Vec2& Vec2::operator*= (const double v1)
{
  Vec2 result{ *this * v1};
  *this = result;
  return *this;
}

Vec2 operator+(const Vec2 &v1, const Vec2 &v2)
{
  return
  {
    v1.x + v2.x,
    v1.y + v2.y
  };
}

Vec2 operator-(const Vec2 &v1, const Vec2 &v2)
{
  return
  {
    v1.x - v2.x,
    v1.y - v2.y
  };
}

Vec2 operator*(const Vec2 &v1, const double scal)
{
  return
  {
    v1.x * scal,
    v1.y * scal
  };
}

Vec2 operator*(const double scal, const Vec2 &v1)
{
  return v1*scal;
}

Vec2 operator/(const Vec2 &v1, const double scal)
{
  if ( scal == 0)
    throw "Division by zero";
  return v1 * (1/scal);
}

Vec2 operator- (const Vec2 &v1)
{
  return
  {
    - v1.x,
    - v1.y
  };
}


Int3 operator+(const Int3 &v1, const Int3 &v2)
{
  return
  {
    v1.x + v2.x,
    v1.y + v2.y,
    v1.z + v2.z
  };
}

Int3& Int3::operator+= (const Int3 &v1)
{
  Int3 result{ *this + v1};
  *this = result;
  return *this;
}

Int3 operator-(const Int3 &v1, const Int3 &v2)
{
  return
  {
    v1.x - v2.x,
    v1.y - v2.y,
    v1.z - v2.z
  };
}

Int3 operator- (const Int3 &v1)
{
  return
  {
    - v1.x,
    - v1.y,
    - v1.z
  };
}

Int3& Int3::operator-= (const Int3 &v1)
{
  return *this+=(-v1);
}

Int3 operator*(const Int3 &v1, const int scal)
{
  return
  {
    v1.x * scal,
    v1.y * scal,
    v1.z * scal
  };
}

Int3 operator*(const int scal, const Int3 &v1)
{
  return v1*scal;
}

Int3& Int3::operator*= (const int v1)
{
  Int3 result{ *this * v1 };
  *this = result;
  return *this;
}

Uint3 operator+(const Uint3 &v1, const Uint3 &v2)
{
  return
  {
    v1.x + v2.x,
    v1.y + v2.y,
    v1.z + v2.z
  };
}

Uint3& Uint3::operator+= (const Uint3 &v1)
{
  Uint3 result{ *this + v1};
  *this = result;
  return *this;
}

Uint3 operator-(const Uint3 &v1, const Uint3 &v2)
{
  return
  {
    v1.x - v2.x,
    v1.y - v2.y,
    v1.z - v2.z
  };
}

Uint3& Uint3::operator-= (const Uint3 &v1)
{
  return *this=*this-v1;
}

Uint3 operator*(const Uint3 &v1, const uint scal)
{
  return
  {
    v1.x * scal,
    v1.y * scal,
    v1.z * scal
  };
}

Uint3 operator*(const uint scal, const Uint3 &v1)
{
  return v1*scal;
}

Uint3& Uint3::operator*= (const uint v1)
{
  Uint3 result{ *this * v1 };
  *this = result;
  return *this;
}

Vec3& Vec3::operator+= (const Vec3 &v1)
{
  Vec3 result{ *this + v1};
  *this = result;
  return *this;
}

Vec3& Vec3::operator-= (const Vec3 &v1)
{
  return *this+=(-v1);
}

Vec3& Vec3::operator*= (const double v1)
{
  Vec3 result{ *this * v1};
  *this = result;
  return *this;
}

Vec3 operator+(const Vec3 &v1, const Vec3 &v2)
{
  return
  {
    v1.x + v2.x,
    v1.y + v2.y,
    v1.z + v2.z
  };
}

Vec3 operator-(const Vec3 &v1, const Vec3 &v2)
{
  return
  {
    v1.x - v2.x,
    v1.y - v2.y,
    v1.z - v2.z
  };
}

Vec3 operator*(const Vec3 &v1, const double scal)
{
  return
  {
    v1.x * scal,
    v1.y * scal,
    v1.z * scal
  };
}

Vec3 operator*(const double scal, const Vec3 &v1)
{
  return v1*scal;
}

Vec3 operator/(const Vec3 &v1, const double scal)
{
  if ( scal == 0)
    throw "Division by zero";
  return v1 * (1/scal);
}

Vec3::~Vec3()
{}

Vec3 operator- (const Vec3 &v1)
{
  return
  {
    - v1.x,
    - v1.y,
    - v1.z
  };
}

std::ostream& operator<< (std::ostream &out, const Vec3 &vec)
{
  out << vec.x << "\t" << vec.y << "\t" << vec.z;
  return out;
}

std::ostream& operator<< (std::ostream &out, const Uint3 &vec)
{
  out << vec.x << "\t" << vec.y << "\t" << vec.z;
  return out;
}

bool operator== (const Vec3 &v1, const Vec3 &v2)
{
  return
  {
    ( v1.x == v2.x ) &&
    ( v1.y == v2.y ) &&
    ( v1.z == v2.z )
  };
}

bool operator!= (const Vec3 &v1, const Vec3 &v2)
{
  return !(v1 == v2);
}

double Vec3::norm() const
{
  return sqrt(dot2(*this));
}

double angularSep(const Vec3 &v1, const Vec3 &v2)
{
  double scalar_product = dot(v1,v2) / ( v1.norm() *v2.norm() );
  double theta{acos(scalar_product)};
  return theta;
}

double getThetaFromOrientation(const Vec3 &orientation)
{
  return atan2(orientation.y,orientation.x);
}

Vec3 min3(Vec3 a, Vec3 b)
{
  return {
    std::min(a.x,b.x),
    std::min(a.y,b.y),
    std::min(a.z,b.z)
  };
}

Vec3 max3(Vec3 a, Vec3 b)
{
  return {
    std::max(a.x,b.x),
    std::max(a.y,b.y),
    std::max(a.z,b.z)
  };
}

Vec3 mod3(Vec3 r, double a)
{
  return {
    std::fmod(r.x,a),
    std::fmod(r.y,a),
    std::fmod(r.z,a)
  };
}

// template signum function
template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

// Return the sign of a double x - This is intended ony for use in this file
inline int sgnDouble(double x)
{
  // Return a value with the magnitude 1 and sign of x
  return sgn(x);
}

double clampRoot(double root)
{
  if (root>=1) return 1;
  else if (root<=-1) return -1;
  else return root;
}

double getPointToLineSegDist(const Vec3 &p, const Vec3 &r1, const Vec3 &r2)
{
  Vec3 v {r2-r1};       // direction vector
  double vv {dot2(v)}; // squared length of line seg

  // if the line is degenerate then just need distance between any end point
  if ( vv==0 ) return (p-r1).norm();

  Vec3 w {p-r1};        // point to line start
  double vw {dot(v,w)}; // projection onto line segment, negative implies obtuse angle

  if (vw < 0)
  {
    // point is before the line start
    return w.norm();
  }
  else if (vv < vw)
  {
    // point is beyond the end of the line
    return (p-r2).norm();
  }
  else
  {
    // the right angle from the point to the line is within the segment
    double mu {vw/vv};  // fraction along seg where the right angle to p intersects
    Vec3 r {r1 + mu * v};
    return (p-r).norm();
  }
}

void closestApproachLineSegmentsPournin(const Vec3 &x_1,
                                        const Vec3 &a_1,
                                        const Vec3 &x_2,
                                        const Vec3 &a_2,
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
  Vec3 x_21 = x_1 - x_2; // vector from x2 to x1
  double a =  dot(a_1,a_1);
  double b = -dot(a_1,a_2);
  double c =  dot(a_2,a_2);
  double d =  dot(a_1,x_21);
  double e = -dot(a_2,x_21);
  double delta = a*c - b*b;  // \equiv modulus( a_1 \cross a_2 )

  // === Handle degenerate cases ===
  if ( a==0 && c==0 ) { s_star=0; t_star=0;               return; }
  if ( a==0 && c>0  ) { s_star=0; t_star=clampRoot(-e/c); return; }
  if ( c==0 && a>0  ) { t_star=0; s_star=clampRoot(-d/a); return; }

  // === Handle line segments ===
  // Handle parallel line segments
  if (delta==0)
  {
    s_star = 0.5*( clampRoot( (c-e)/b )+clampRoot( -(c+e)/b ) );
    t_star = clampRoot(-(e+b*s_star)/c);
    // std::cout << "s: " << s_star << " t: " << t_star << '\n';
  }
  else
  {
    t_star = clampRoot( ( b*d - a*e ) / delta );

    double s_star_tmp = (-b*t_star - d) / a;


    s_star = clampRoot( s_star_tmp );

    if ( s_star_tmp < -1 || s_star_tmp > 1 )
      t_star = clampRoot( ( -b*s_star - e ) / c );
  }
}
