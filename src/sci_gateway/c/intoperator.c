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

//#include <stack-c.h>
#include <api_scilab.h>

#include "intoperator.h"

int
intoperator (char *fname)
{
  int ainf, asup, binf, bsup, op;
  int sa[2], sb[2], flag[1];
  int cinf, csup, f;
  int ma, na, ma2, na2;
  int mb, nb, mb2, nb2;
  int r, s;
  int p, q, t, u;
  int rt;

  CheckRhs (5, 5);
  CheckLhs (2, 3);

  GetRhsVar (1, "d", &na, &ma, &ainf);
  GetRhsVar (2, "d", &na2, &ma2, &asup);
  GetRhsVar (3, "c", &t, &u, &op);
  GetRhsVar (4, "d", &nb, &mb, &binf);
  GetRhsVar (5, "d", &nb2, &mb2, &bsup);

  if (*cstk (op) == 'k' || *cstk (op) == 'y')
    {
      r = na * nb;
      s = ma * mb;
    }

  else
    {
      if (na * ma == 1)
	{
	  r = nb;
	  s = mb;
	}

      else if (nb * mb == 1 || (na == nb && ma == mb))
	{
	  r = na;
	  s = ma;
	}

      else if (ma == nb)
	{
	  r = na;
	  s = mb;
	}

      else
	{
	  r = 0;
	  s = 0;
	}
    }

  p = 1;
  q = 1;

  CreateVar (6, "d", &r, &s, &cinf);
  CreateVar (7, "d", &r, &s, &csup);
  CreateVar (8, "i", &p, &q, &f);

  sa[0] = na;
  sa[1] = ma;
  sb[0] = nb;
  sb[1] = mb;

  if (r * s != 0)
    {
      if (nb2 * mb2 == 0)
	{
	  operate (stk (ainf), stk (asup), sa,
		   stk (binf), NULL, sb,
		   *cstk (op), stk (cinf), stk (csup), flag);
	}

      else if (na2 * ma2 == 0)
	{
	  operate (stk (ainf), NULL, sa,
		   stk (binf), stk (bsup), sb,
		   *cstk (op), stk (cinf), stk (csup), flag);
	}

      else
	{
	  operate (stk (ainf), stk (asup), sa,
		   stk (binf), stk (bsup), sb,
		   *cstk (op), stk (cinf), stk (csup), flag);
	}
    }

  else
    {
      flag[0] = -1;
    }

  *istk (f) = flag[0];

  LhsVar (1) = 6;
  LhsVar (2) = 7;
  LhsVar (3) = 8;

  return 0;

}

int
intpower (char *fname)
{
  int ainf, asup, op, n;
  int sa[2], flag[1];
  int cinf, csup, f;
  int ma, na;
  int r, s, v, w;
  int p, q, t, u;

  CheckRhs (4, 4);
  CheckLhs (2, 3);

  GetRhsVar (1, "d", &na, &ma, &ainf);
  GetRhsVar (2, "d", &na, &ma, &asup);
  GetRhsVar (3, "c", &t, &u, &op);
  GetRhsVar (4, "i", &p, &q, &n);

  r = na;
  s = ma;
  v = 1;
  w = 1;

  CreateVar (5, "d", &r, &s, &cinf);
  CreateVar (6, "d", &r, &s, &csup);
  CreateVar (7, "i", &v, &w, &f);

  sa[0] = na;
  sa[1] = ma;

  power (stk (ainf), stk (asup), sa, *cstk (op), istk (n),
	 stk (cinf), stk (csup), flag);

  *istk (f) = flag[0];

  LhsVar (1) = 5;
  LhsVar (2) = 6;
  LhsVar (3) = 7;

  return 0;

}

int
intunary (char *fname)
{
  int ainf, asup, op;
  int sa[2], flag[1];
  int cinf, csup, f;
  int ma, na;
  int r, s;
  int t, u, v, w;

  CheckRhs (3, 3);
  CheckLhs (2, 3);

  GetRhsVar (1, "d", &na, &ma, &ainf);
  GetRhsVar (2, "d", &na, &ma, &asup);
  GetRhsVar (3, "c", &t, &u, &op);

  r = na;
  s = ma;
  v = 1;
  w = 1;

  CreateVar (4, "d", &r, &s, &cinf);
  CreateVar (5, "d", &r, &s, &csup);
  CreateVar (6, "i", &v, &w, &f);

  sa[0] = na;
  sa[1] = ma;

  unary (stk (ainf), stk (asup), sa, *cstk (op),
	 stk (cinf), stk (csup), flag);

  *istk (f) = flag[0];

  LhsVar (1) = 4;
  LhsVar (2) = 5;
  LhsVar (3) = 6;

  return 0;

}

int
introundmode (char *fname)
{
  int direct;
  int t, u, v;

  CheckRhs (1, 1);
  CheckLhs (1, 1);
  
  GetRhsVar (1, "i", &t, &u, &direct);
  CreateVar (2, "i", &t, &u, &v);
  
  *istk (v) = roundmode(istk (direct)); 
  
  LhsVar (1) = 2;
  
  return 0;
}
