/*----------------------------------------------------------------------------*/
/*
  MIT License

  Copyright (c) 2019 Mya Warren, Hui Sun, Yue Yan, Bo Li, Terry Hwa

  Permission is hereby granted, free of charge, to any person obtaining a copy
  of this software and associated documentation files (the "Software"), to deal
  in the Software without restriction, including without limitation the rights
  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
  copies of the Software, and to permit persons to whom the Software is
  furnished to do so, subject to the following conditions:

  The above copyright notice and this permission notice shall be included in all
  copies or substantial portions of the Software.

  Modification from the original by Rory Claydon 2021

*/
#ifndef GEOMETRY_HPP
#define GEOMETRY_HPP

// Custom modules
#include "IBacterium.hpp"
#include "MathUtility.hpp"        // definition of Vec3 class

// Enforce a periodic boundary condition for a box of size L in the domain
// [-L/2,L/2]
// Extra L is to catch particles moving more negative than they should
inline double periodicBC(double x, double L)
{
	return std::fmod(x+1.5*L,L) - 0.5*L;
}
inline Vec3 periodicBC(Vec3 x, double L)
{
	return mod3(x+1.5*L,L) - 0.5*L;
}

inline double clamp(double n, double minn, double maxn)
{
	if (n<minn)
		return minn;
	else if (n>maxn)
		return maxn;
	else
		return n;
}

inline double clampRoot(double root)
{
  if (root>=1) return 1;
  else if (root<=-1) return -1;
  else return root;
}

inline bool checkContact(
	IBacterium* cellA,
	IBacterium* cellB,
	double &s,
	double &t,
	Vec3 &cv
)
{
	const Vec3 v1 { 0.5*cellA->getOrientation()*cellA->getLength() };
	const Vec3 v2 { 0.5*cellB->getOrientation()*cellB->getLength() };
	closestApproachLineSegmentsPournin(cellA->getPos(),v1,
																		 cellB->getPos(),v2,
																		 s,t);
	const Vec3 ca { cellA->getPos()+v1*s };
	const Vec3 cb { cellB->getPos()+v2*t };
  cv = cb-ca;
  return dot2(cv) <= dot2(cellA->getRadius()+cellB->getRadius());
}


/**
  @details
  Minimise the distance between two line segments.
  Define R(s,t) = a s^2 + 2 b s t + c t^2 + 2 d s + 2 e t + f
                = |x_1 - x_2|^2
  Return the values s^* and t^* s.t. R(s^*,t^*) is minimal subject to the
  constraint that s,t \in [-1,1]^2

  If the lines are parallel then the proximate points are set to be the
  midpoint of the overlap.

  Algorithm from [1] and [2]
*/
inline bool closestApproachLineSegmentsParallel(
  const Vec3 x_1,
  const Vec3 a_1,
  const Vec3 x_2,
  const Vec3 a_2,
  double &s1,
  double &t1,
  double &s2,
  double &t2
)
{
	constexpr double EPSILON_SQR { 1e-4 };
  Vec3 x_21 = x_1 - x_2; // vector from x2 to x1
  double a =  dot(a_1,a_1);
  double b = -dot(a_1,a_2);
  double c =  dot(a_2,a_2);
  double d =  dot(a_1,x_21);
  double e = -dot(a_2,x_21);
  double delta = a*c - b*b; //!< \sim modulus( \hat{a}_1 \cross \hat{a}_2 )^2

#ifdef SPHERICAL // Avoid branches if spherical particles are not present
  // === Handle degenerate cases ===
  if ( a==0 && c==0 ) { s1=0; t1=0;               return; }
  if ( a==0 && c>0  ) { s1=0; t1=clampRoot(-e/c); return; }
  if ( c==0 && a>0  ) { t1=0; s1=clampRoot(-d/a); return; }
#endif

  // === Handle line segments ===
  if ( delta<EPSILON_SQR*a*c )       // Segments deemed close enough to parallel
  {
    s1 = clampRoot( (c-e)/b );
    t1 = clampRoot(-(e+b*s1)/c);
    s2 = clampRoot( -(c+e)/b );
    t2 = clampRoot(-(e+b*s2)/c);
    return true;
  }
  else
  {
    t1 = clampRoot( ( b*d - a*e ) / delta );

    double s1_tmp = (-b*t1 - d) / a;

    s1 = clampRoot( s1_tmp );

    if ( s1_tmp < -1 || s1_tmp > 1 )
      t1 = clampRoot( ( -b*s1 - e ) / c );
    return false;
  }
}

inline void getMinDist(
	const IBacterium* cell1,
  const IBacterium* cell2,
  double& d,
  Vec3& c1,
  Vec3& c2
)
{
	// Line segments end points
	Vec3 p1,q1,p2,q2;
	cell1->getMyEndVecs(p1,q1);
	cell2->getMyEndVecs(p2,q2);

	// find the direction Segment of S1 and S2
	const Vec3 v1 = 0.5*(q1-p1);
	const Vec3 v2 = 0.5*(q2-p2);

	// s,t parameterise the line segments
	double s,t;
	closestApproachLineSegmentsPournin(cell1->getPos(),v1,
																		 cell2->getPos(),v2,
																		 s, t);
	c1 = cell1->getPos()+v1*s;
	c2 = cell2->getPos()+v2*t;
	const Vec3 cv = c1-c2;
	d = sqrt(dot(cv,cv));		// length of segment connecting two cells
}

// find the minimum distance between any two cells
inline void getMinDistWarren(
  const IBacterium* cell1,
  const IBacterium* cell2,
  double& d,
  Vec3& c1,
  Vec3& c2
)
{
	Vec3 p1,q1,p2,q2,v1,v2,r,cv;
	double l1,l2,f,c,b,denom,tnom,s,t;

  // S1 and S2 are line Segments defined by the points (p1,q1) and (p2,q2)
  // respectively
  cell1->getMyEndVecs(p1,q1);
	cell2->getMyEndVecs(p2,q2);

	// find the direction Segment of S1 and S2
	v1 = q1-p1;
	v2 = q2-p2;
	r  = p1-p2;	// vec between two start points from 2 to 1

	l1 = dot(v1,v1); // length of segment 1
	l2 = dot(v2,v2); // length segment 2
	f  = dot(v2,r);  // projection of segment 2 in the direction of 1
	c  = dot(v1,r);  // projection of segment 1 in the same vain
	b  = dot(v1,v2); // projection of v1 along v2

	if ( l1==0 && l2==0 ) // both segments are degenerate
	{
		s=t=0;
	}
	else if ( l1==0 ) // first segment is just a point
	{
		s=0;
		t = clamp(f/l2,0,1);
	}
	else if ( l2==0 ) // second segment is just a point
	{
		s=clamp(-c/l1,0,1);
		t=0;
	}
	else // both segments are not degenerate
	{
		denom = l1*l2-b*b; // | v1 x v2 |^2

		if (denom!=0) // if not parallel Segments or degenerate
		// closest point on cell1 to cell2 within the Segment cell1
		// try to find the vector between the infinite line segments then clamp it
		{
			s = clamp((b*f-c*l2)/denom,0,1);
		}
		else
		{
			// If parallel select the midpoint of the overlap (coming soon!)
			s = 0;
			std::cout << "\n\n----------- parallel -----------\n\n";
		}

		// compute point on L2 closest to S1

		tnom = b*s+f; // first guess for t

	  // check if this is within the Segment (t between 0 and 1)
		if (tnom<0)
		{
			t = 0;
			assert(l1>0);
			s = clamp(-c/l1,0,1);
		}
		else if (tnom>l2)
		{
		  t = 1;
			assert(l1>0);
		  s = clamp((b-c)/l1,0,1);
		}
		else
		{
			assert(l2>0);
		  t = tnom/l2;
		}
	}
	c1 = p1+v1*s;
	c2 = p2+v2*t;
	cv = c1-c2;
	d = sqrt(dot(cv,cv));		// length of segment connecting two cells
}

inline void findSqrDistToLine(
  Vec3  r,
  Vec3  n,
  float d,
  Vec3  &rr,
  double &dist2,
  double &t
)
{
  Vec3 dr=rr-r;
  t=dot(dr,n);
  if (t>d) t=d ;
  if (t<-d) t=-d ;
  dist2=dot2(dr-t*n);
}

inline void closestApproachLineSegments(
  const Vec3 ri,
  Vec3 ni,
  const double  di,
  const Vec3 rj,
  Vec3 nj,
  const double  dj,
  double &ti,
  double &tj,
  double &dij
)
{
  // Find which end is closest
  Vec3 rij { ri-rj };
  double ninj=dot(ni,nj), w=ninj*ninj-1;
  // w is -|ni x nj|^2
  if (w==0) { ti=tj=0; }
  else {
    double snirij=dot(ni,rij),
    snjrij=dot(nj,rij);
    ti=( snirij-ninj*snjrij)/w;
    tj=(-snjrij+ninj*snirij)/w;
    if (ti>di)  ti= di;
    if (ti<-di) ti=-di;
    if (tj>dj)  tj= dj;
    if (tj<-dj) tj=-dj;
  }
  ni*=ti; nj*=tj;
  Vec3 d=rij; d+=ni; d-=nj;
  dij=dot2(d);
}

template <class P>
inline void getMinDistBartek(
  const P *A,
  const P *B,
	double &dmin,
  Vec3 &cA,
  Vec3 &cB
)
{
  Vec3 r1,r2,ri1,ri2,rj1,rj2,bin,bjn,dr1,dr2,pos_A, pos_B;
  double di1, di2, dj1, dj2, ti1, ti2, tj1, tj2, dij, ti, tj;
  double di { 0.5*A->getLength() };
  double dj { 0.5*B->getLength() };
  pos_A=A->getPos();
  pos_B=B->getPos();
  Vec3 n_A = A->getOrientation();
  Vec3 n_B = B->getOrientation();
  closestApproachLineSegments(
    pos_A,
    n_A,
    di,
    pos_B,
    n_B,
    dj,
    ti,
    tj,
    dij
  );
  bin=n_A; bin*=di;
  bjn=n_B; bjn*=dj;
  ri1=pos_A; ri1+=bin;
  ri2=pos_A; ri2-=bin;
  rj1=pos_B; rj1+=bjn;
  rj2=pos_B; rj2-=bjn;

  findSqrDistToLine(pos_B,n_B,dj,ri1,di1,ti1);
  findSqrDistToLine(pos_B,n_B,dj,ri2,di2,ti2);
  findSqrDistToLine(pos_A,n_A,di,rj1,dj1,tj1);
  findSqrDistToLine(pos_A,n_A,di,rj2,dj2,tj2);
  uint flag{ 0 };
  if      (di1<=di2 && di1<=dj1 && di1<=dj2 && di1<=dij) { dmin=di1; flag=1; r1=ri1; r2=n_B; r2*=ti1; r2+=pos_B; }
  else if (di2<=di1 && di2<=dj1 && di2<=dj2 && di2<=dij) { dmin=di2; flag=1; r1=ri2; r2=n_B; r2*=ti2; r2+=pos_B; }
  else if (dj1<=di1 && dj1<=di2 && dj1<=dj2 && dj1<=dij) { dmin=dj1; flag=2; r1=rj1; r2=n_A; r2*=tj1; r2+=pos_A; }
  else if (dj2<=di1 && dj2<=di2 && dj2<=dj1 && dj2<=dij) { dmin=dj2; flag=2; r1=rj2; r2=n_A; r2*=tj2; r2+=pos_A; }
  else if (dij<=di1 && dij<=di2 && dij<=dj1 && dij<=dj2) { dmin=dij; flag=1; r1=pos_A; r1+=n_A*ti; r2=pos_B; r2+=n_B*tj; }
  if ( flag==1 ) {  dr1=r1; dr1-=pos_A; dr2=r2; dr2-=pos_B; cA=r1; cB=r2; }
  if ( flag==2 ) {  dr1=r1; dr1-=pos_B; dr2=r2; dr2-=pos_A; cA=r2; cB=r1; }
	dmin=sqrt(dmin);
}
#endif // end fileguard
