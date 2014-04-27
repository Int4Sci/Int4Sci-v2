//    <Int4Sci, a scilab interface for Interval Analysis.>
//    Copyright (C) <2007>  <R. PEREIRA, D. DANEY, J.P. MERLET>

//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.

//   This program is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.

//   You should have received a copy of the GNU General Public License
//   along with this program.  If not, see <http://www.gnu.org/licenses/>.


//   Contact : R. PEREIRA <int4sci@lists-sop.inria.fr>
//   COPRIN team
//   INRIA Sophia Antipolis
//   2004 route des lucioles - BP 93
//   06902 Sophia Antipolis Cedex
//   FRANCE

#include <stack-c.h>

#include "intILSR.h"

int
intILSR (char *fname)
{

  int ainf, asup, binf, bsup, xinf, xsup, op;
  int sa[2], sb[2], sx[2], flag[1];
  int cinf, csup, f;
  int ma, na;
  int mb, nb;
  int mx, nx;
  int r, s, o, v;
  int p, q, t, u;

  CheckRhs (5, 7);
  CheckLhs (2, 3);

  GetRhsVar (1, "d", &na, &ma, &ainf);
  GetRhsVar (2, "d", &na, &ma, &asup);
 
  GetRhsVar (3, "d", &nb, &mb, &binf);
  GetRhsVar (4, "d", &nb, &mb, &bsup);

  GetRhsVar (5, "d", &nx, &mx, &xinf);
  GetRhsVar (6, "d", &nx, &mx, &xsup);

  GetRhsVar (7, "c", &t, &u, &op);

  r = nb;
  s = mb;

  o = nx;
  v = mx;

  p = 1;
  q = 1;

  CreateVar (8, "d", &r, &s, &cinf);
  CreateVar (9, "d", &r, &s, &csup);
  CreateVar (10, "i", &p, &q, &f);

  sa[0] = na;
  sa[1] = ma;
  sb[0] = nb;
  sb[1] = mb;
 
  if (o * v != 0)
    {
      sx[0] = nx;
      sx[1] = mx;

      ILSR (stk (ainf), stk(asup), sa,
	    stk (binf), stk(bsup), sb,
	    stk (xinf), stk(xsup), sx,
	    *cstk (op), stk (cinf), stk (csup), flag);

    }

  else
    {
      ILSR (stk (ainf), stk(asup), sa,
	    stk (binf), stk(bsup), sb,
	    NULL, NULL, sb,
	    *cstk (op), stk (cinf), stk (csup), flag);
    }

  *istk (f) = flag[0];

  LhsVar (1) = 8;
  LhsVar (2) = 9;
  LhsVar (3) = 10;

  return 0;

}

int 
intI4Svarsend ( char *fname )
{
  int a, b, op;
  int ma, na;
  int mb, nb;
  int t, u;
  int v, w, x, y;

  CheckRhs (1, 3);
  CheckLhs (0, 1);

  GetRhsVar (1, "c", &t, &u, &op);

  GetRhsVar (2, "d", &na, &ma, &a);
  GetRhsVar (3, "i", &nb, &mb, &b);

  v = na;
  w = ma;

  x = nb;
  y = mb;

  if ( x*y == 0 )
  {
	I4Sreceived(*cstk(op), stk(a), NULL);
  }
  
  else if ( v*w == 0 )
  {
	I4Sreceived(*cstk(op), NULL, istk(b));
  }

  else
  {
	I4Sreceived(*cstk(op), stk(a), istk(b));	
  }

  return 0;

}


int 
intI4Svarget ( char *fname )
{
	int op, c, d;
	int t, u, r, s;
	int m, n;

	CheckRhs (1, 1);
	CheckLhs (1,2);

	GetRhsVar (1, "c", &t, &u, &op);
	
	r = 1;
	s = 1;
	m = 1;
	n = 1;

	CreateVar (2, "d", &r, &s, &c);
	CreateVar (3, "i", &m, &n, &d);

	I4Ssent(*cstk (op), stk (c), istk (d) );

	LhsVar (1) = 2;
	LhsVar (2) = 3;

}

