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


/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 * HELP
 * This file is divided in 4 sub part
 * - I the general functions interfaced in scilab
 * - II the preconditionning
 * - III Krawcyck, Gauss elimination, Gauss sidel, Bareiss ,...
 * - IV Tools used ..

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#include "ProfilILSR.h"
#include "intILSR.h"
using namespace std;

int MAX_LOOP;
double FAC_IMPROVE;

#define iterV2M(k,nc) (k/nc)+1,(k%nc)+1
#define iterf2c(k,nl,nc) ((k%nc)*nl+(k/nc))
#define iterc2f(k,nl,nc) ((k%nl)*nc+(k/nl))



/*######################
 * I General functions
 *#####################*/

extern "C"
{
  int
    ILSR (const double *ainf, const double *asup, const int *sa,
	  const double *binf, const double *bsup, const int *sb,
	  const double *xinf, const double *xsup, const int *sx,
	  char op, double *cinf, double *csup, int *flag)
  {
   /*---------------------------------------
   * Entry flag conventions for resolution method choice
   * 'a' : Profil LSS (Not implemented)
   * 'b' : PRofil ILSS
   * 'c' : classical Krawczyck 
   * 'd' : high level Krawczyck
   * 'e' : Gauss Elimination
   * 'f' : Bareiss
   * 'g' : Gauss Sidel
   * 'h' : LU
   * 'i' : Hansen Bliek
   * 'j' : classical Krawczyck with preconditionning
   * 'k' : high level Krawczyck with preconditionning
   * 'l' : Gauss Elimination with preconditionning
   * 'm' : Bareiss with preconditionning
   * 'n' : Gauss Sidel with preconditionning
   * 'o' : LU with preconditionning
   * 'p' : Hansen Bliek with preconditionning
   *
   * Exit flag conventions
   * >= -3 : Not Implented Yet or precondition problems
   * -2 : not well constraint case
   * -1 or 0 : there's a problem on the tests made during
   *            the method execution or no solution
   *  1 : finish without problem
   ----------------------------------------*/

    INTERVAL_MATRIX pA (sa[0], sa[1]);
    INTERVAL_VECTOR pB (sb[1] * sb[0]);
    INTERVAL_VECTOR X (sx[1] * sx[0]);
    VECTOR CI (sb[1] * sb[0]);
    VECTOR CS (sb[1] * sb[0]);
    int OK;
    int i;

    for (i = 0; i < sa[0] * sa[1]; i++)
      {
	INTERVAL A (ainf[iterf2c (i, sa[0], sa[1])],
		    asup[iterf2c (i, sa[0], sa[1])]);
	  pA (iterV2M (i, sa[1])) = A;
      }

    for (i = 0; i < sb[0] * sb[1]; i++)
      {
	INTERVAL B (binf[iterf2c (i, sb[0], sb[1])],
		    bsup[iterf2c (i, sb[0], sb[1])]);
	pB (i + 1) = B;
      }

    if (xinf != NULL)
      {
	for (i = 0; i < sx[0] * sx[1]; i++)
	  {
	    INTERVAL Y (xinf[iterf2c (i, sx[0], sx[1])],
			xsup[iterf2c (i, sx[0], sx[1])]);
	    X (i + 1) = Y;
	  }
      }

    /* --------------------------------------
       Profil ILSS : I4S default execution
       ----------------------------------------- */

    *flag = -3;

    if (op == 'b')
      {
	X = ILSS (pA, pB, OK);

	if (OK == 1)
	  {
	    *flag = 1;
	  }

	else
	  {
	    *flag = -1;
	  }
      }

  /*--------------------------
    No initial condition on X 
  ----------------------------*/

    if (op == 'e')
      {
	*flag = ilsge (pA, pB, X);
      }

    if (op == 'f')
      {
	*flag = ilsb (pA, pB, X);
      }

    if (op == 'h')
      {
	*flag = -3;
      }

    if (op == 'i')
      {
	*flag = ilshb (pA, pB, X);
      }

    if (op == 'l')
      {
	*flag = pils (ilsge, pA, pB, X);
      }

    if (op == 'm')
      {
	*flag = pils (ilsb, pA, pB, X);
      }

    if (op == 'o')
      {
	*flag = -3;
      }

    if (op == 'p')
      {
	*flag = pils (ilshb, pA, pB, X);
      }

  /*---------------------
    Initial condition on X 
  --------------------------*/

    if (op == 'c')
      {
	*flag = ilsk (pA, pB, X);
      }

    if (op == 'd')
      {
	*flag = ilskgs (pA, pB, X);
      }

    if (op == 'g')
      {
	*flag = ilsgs (pA, pB, X);
      }

    if (op == 'j')
      {
	*flag = pils (ilsk, pA, pB, X);
      }

    if (op == 'k')
      {
	*flag = pils (ilskgs, pA, pB, X);
      }

    if (op == 'n')
      {
	*flag = pils (ilsgs, pA, pB, X);
      }

    CI = Inf (X);
    CS = Sup (X);

    for (i = 0; i < sb[0] * sb[1]; i++)
      {
	cinf[iterf2c (i, sb[0], sb[1])] = CI (i + 1);
	csup[iterf2c (i, sb[0], sb[1])] = CS (i + 1);
      }

    return *flag;
  }

}

/*###############################
Changing initial constant values
##################################*/

int
I4Sreceived (char var, double *fval, int *val)
{
  if (var == '0')
    {
      FAC_IMPROVE = fval[0];
      MAX_LOOP = val[0];
    }

  if (var == '1')
    {
      MAX_LOOP = val[0];
    }

  if (var == '2')
    {
      FAC_IMPROVE = fval[0];
    }

  return 1;
}

int
I4Ssent (char var, double *fval, int *val)
{
  if (var == '0')
    {
      fval[0] = FAC_IMPROVE;
      val[0] = MAX_LOOP;
    }

  if (var == '1')
    {
      val[0] = MAX_LOOP;
    }

  if (var == '2')
    {
      fval[0] = FAC_IMPROVE;
    }

  return 1;
}

/*######################
 * II Preconditionner
 *#####################*/

int
pils (fctILS prog, const INTERVAL_MATRIX & A, const INTERVAL_VECTOR & b,
      INTERVAL_VECTOR & x)
{
  /*---------------------------------------
   * Preconditionning Algorithm by middle of A (c=[mid(A)]^-1): 
   * give x such that C.A.x=C.b
   * Algorithm described in 
   * ???
   * RETURN
   * -3 under-constraint case
   * -4 the square mid point matrix is singular 
   * -5 the square mid point matrix is singular 
   * else return the ouput of prog
   ----------------------------------------*/
  INT nbrow = RowDimension (A);
  INT nbcol = ColDimension (A);
  MATRIX C (nbcol, nbrow);
  /*****************************************
   * under-constraint problem 
   ******************************************/
  if (nbrow < nbcol)
    return -3;
  /*****************************************
   * well-constraint problem 
   ******************************************/
  if (nbrow = nbcol)
    {
      if (MatInverse (Mid (A), C) == -1)
	if (MatInverse (Inf (A), C) == -1)
	  if (MatInverse (Sup (A), C) == -1)
	    return -4;
      return prog (C * A, C * b, x);
    }

  /*****************************************
   * over-constraint problem 
   ******************************************/
  if (nbrow > nbcol)
    {
      MATRIX Ct (nbcol, nbcol);
      MATRIX mA (nbrow, nbcol);
      MATRIX mAt (nbcol, nbrow);
      mA = Mid (A);
      mAt = Transpose (mA);
      if (MatInverse (mAt * mA, Ct) == -1)
	return -5;		/* pseudo-inverse */
      C = Ct * mAt;
      return prog (C * A, C * b, x);
    }
}

/*################################################
 * III Interval linear solving classical algorithm
 *################################################*/

/*------------------------------------
 *             Krawczyck classique
 *-----------------------------------*/

int
ilsk (const INTERVAL_MATRIX & A, const INTERVAL_VECTOR & b,
      INTERVAL_VECTOR & x)
{
  /*---------------------------------------
   * Krawczyck Algorithm : give x such that A.x=b
   * Algorithm described in 
   * "Introduction to numerical Analysis"
   * A. Neumaier, cambridge university press 2001
   * RETURN
   * -2 not well constraint case
   * -1 Krawczyck's test faile
   * 0 non intersection during iterative Krawczyk schem 
   * 1 finish without problem
   ----------------------------------------*/
  INT nbrow = RowDimension (A);
  INT nbcol = ColDimension (A);
  /*****************************************
   * not well-constraint problem 
   ******************************************/
  if (nbrow != nbcol)
    return -2;
  /***********************************
   * Temporary Matrix and Vector need
   ************************************/
  INTERVAL_MATRIX K (nbcol, nbcol);
  K = Id (nbcol) - A;
  /********************
   * Krawczyck's test
   ********************/
  REAL kt = INorm (K);
  if (kt > 1.0)
    return -1;
  /****************************
   * Algorithm
   ****************************/
  int iter = 1;
  REAL vol_old = Machine::PosInfinity;
  REAL vol = Perimeter (x);
  REAL fac = (1 + kt) / 2;
  INTERVAL_VECTOR xt (Dimension (x));
  Clear (xt);
  if (fac > FAC_IMPROVE)
    fac = FAC_IMPROVE;
  while (vol <= fac * vol_old && iter++ <= MAX_LOOP)
    {
      xt = (b + K * x);
      if (Intersection (xt, x, xt) == 0)
	return 0;
      x = xt;
      vol_old = vol;
      vol = Perimeter (x);
    }
  return 1;
}

/*------------------------------------
 *   Krawczyck 
 *  direct actuallization of unknowns 
 *-----------------------------------*/

int
ilskgs (const INTERVAL_MATRIX & A, const INTERVAL_VECTOR & b,
	INTERVAL_VECTOR & x)
{
  /*---------------------------------------
   * Krawczyck Algorithm : give x such that A.x=b
   * Algorithm described in 
   * "Introduction to numerical Analysis"
   * A. Neumaier, cambridge university press 2001
   * RETURN
   * -2 not well-constraint case
   * -1 Krawczyck's test fail
   * 0 non intersection during iterative Krawczyk schem 
   * 1 finish without problem
   ----------------------------------------*/
  INT nbrow = RowDimension (A);
  INT nbcol = ColDimension (A);
  /*****************************************
   * not well-constraint problem 
   ******************************************/
  if (nbrow != nbcol)
    return -2;
  /***********************************
   * Temporary Matrix and Vector need
   ************************************/
  INTERVAL_MATRIX K (nbrow, nbcol);
  K = Id (nbcol) - A;
  /********************
   * Krawczyck's test
   ********************/
  REAL kt = INorm (K);
  if (kt >= 1.0)
    return -1;
  /****************************
   * Algorithm
   ****************************/
  int iter = 1;
  REAL vol_old = Machine::PosInfinity;
  REAL vol = Perimeter (x);
  REAL fac = (1 + kt) / 2;
  INTERVAL tmp, xt;
  if (fac > FAC_IMPROVE)
    fac = FAC_IMPROVE;
  while (vol <= fac * vol_old && iter++ <= MAX_LOOP)
    {
      for (int i = 1; i <= nbrow; i++)
	{
	  tmp = Hull (0);
	  for (int j = 1; j <= nbcol; j++)
	    tmp += K (i, j) * x (j);
	  xt = b (i) + tmp;
	  if (Intersection (xt, x (i), xt) == 0)
	    return 0;
	  x (i) = xt;
	}
      vol_old = vol;
      vol = Perimeter (x);
    }
  return 1;
}

/*------------------------------------
 *       Gauss Elimination
 *-----------------------------------*/

int
ilsge (const INTERVAL_MATRIX & Ain, const INTERVAL_VECTOR & bin,
       INTERVAL_VECTOR & x)
{
  /*---------------------------------------
   * Gauss Elimination Algorithm : give x such that A.x=b
   * Algorithm described in 
   * ""
   * L. Jaulin and al
   * RETURN
   * -2 not well-constraint problem 
   * -1 null pivot
   * 0 
   * 1 finish without problem
   ----------------------------------------*/

  if (RowDimension (Ain) != ColDimension (Ain))
    return -2;
  int dim = RowDimension (Ain);
  INTERVAL_MATRIX A (dim, dim);
  INTERVAL_VECTOR b (dim);
  INTERVAL_VECTOR tmp (dim);
  INTERVAL sum;
  int i, j, k;
  A = Ain;
  b = bin;
  for (i = 1; i <= dim - 1; i++)
    {
      if (Sup (A (i, i)) >= 0 && Inf (A (i, i)) <= 0)
	return -1;
      for (j = i + 1; j <= dim; j++)
	{
	  tmp (j) = A (j, i) / A (i, i);
	  b (j) = b (j) - tmp (j) * b (i);
	  for (k = i + 1; k <= dim; k++)
	    A (j, k) = A (j, k) - tmp (j) * A (i, k);
	}
    }
  for (i = dim; i >= 1; i--)
    {
      sum = b (i);
      for (j = i + 1; j <= dim; j++)
	sum -= A (i, j) * b (j);
      if (Sup (A (i, i)) >= 0 && Inf (A (i, i)) <= 0)
	return -1;
      b (i) = sum / A (i, i);
    }
  x = b;
  return 1;
}

/*------------------------------------
  A=Ain;
  b=bin;
  for(i=1;i<=dim-1;i++){
    if(Sup(A(i,i))>=0 && Inf(A(i,i))<=0) return -1;
    for(j=i+1;j<=dim;j++){
      b(j)=b(j)*A(i,i)-A(j,i)*b(i);
      for(k=i+1;k<=dim;k++)

 *       Bareiss
 *-----------------------------------*/

int
ilsb (const INTERVAL_MATRIX & Ain, const INTERVAL_VECTOR & bin,
      INTERVAL_VECTOR & x)
{
  /*---------------------------------------
   * Bareis Algorithm : give x such that A.x=b
   * Algorithm described in 
   * "??"
   * 
   * RETURN
   * -2 not well-constraint problem 
   * -1 null pivot
   * 0 
   * 1 finish without problem
   ----------------------------------------*/

  if (RowDimension (Ain) != ColDimension (Ain))
    return -2;

  int dim = RowDimension (Ain);
  INTERVAL_MATRIX A (dim, dim);
  INTERVAL_VECTOR b (dim);
  INTERVAL_VECTOR tmp (dim);
  INTERVAL sum;
  int i, j, k;
  A = Ain;
  b = bin;
  for (i = 1; i <= dim - 1; i++)
    {
      if (Sup (A (i, i)) >= 0 && Inf (A (i, i)) <= 0)
	return -1;
      for (j = i + 1; j <= dim; j++)
	{
	  b (j) = b (j) * A (i, i) - A (j, i) * b (i);
	  for (k = i + 1; k <= dim; k++)
	    A (j, k) = A (j, k) * A (i, i) - A (j, i) * A (i, k);
	}
    }
  for (i = dim; i >= 1; i--)
    {
      sum = b (i);
      for (j = i + 1; j <= dim; j++)
	sum -= A (i, j) * b (j);
      if (Sup (A (i, i)) >= 0 && Inf (A (i, i)) <= 0)
	return -1;
      b (i) = sum / A (i, i);
    }
  x = b;
  return 1;
}

/*------------------------------------
 *       Gauss Sidel
 *-----------------------------------*/

int
ilsgs (const INTERVAL_MATRIX & A, const INTERVAL_VECTOR & b,
       INTERVAL_VECTOR & x)
{
  /*---------------------------------------
   * Gauss Sidel Algorithm : give x such that A.x=b
   * Algorithm described in 
   * ""
   * L. Jaulin and al
   * RETURN
   * -2 not well-constraint problem 
   * -1 No solution
   * 0 
   * 1 finish without problem
   ----------------------------------------*/
  INT nbrow = RowDimension (A);
  INT nbcol = ColDimension (A);
  if (nbrow != nbcol)
    return -2;
  int dim = nbcol;
  int iter = 1;
  REAL vol_old = Machine::PosInfinity;
  REAL vol = Perimeter (x);
  INTERVAL tmp, temp;
  int i, j, k;
  INTERVAL before, after, xt;

  while (vol <= FAC_IMPROVE * vol_old && iter++ <= MAX_LOOP)
    {
      for (i = 1; i <= dim; i++)
	{
	  if (Sup (A (i, i)) < 0 || Inf (A (i, i)) > 0)
	    {
	      before = Hull (0.0);
	      after = Hull (0.0);
	      for (j = 1; j < i; j++)
		{
		  before += A (i, j) * x (j);
		}
	      for (j = i + 1; j <= dim; j++)
		{
		  after += A (i, j) * x (j);
		}
	      if (Intersection
		  (xt, x (i), (b (i) - after - before) / A (i, i)) == 0)
		return -1;
	      if (Inf (xt) == Inf (x (i)) && Sup (xt) == Sup (x (i)))
		return -1;
	      x (i) = xt;
	    }
	  vol_old = vol;
	  vol = Perimeter (x);
	}
    }
  return 1;
}

int
ilshb (const INTERVAL_MATRIX & A, const INTERVAL_VECTOR & B,
       INTERVAL_VECTOR & x)
{
/*--------------------------------------------------
 *Hansen-Bliek Algorithm : 
 * Compute hull of SIGMA(A,b) if A is identity-centered and inf(A) is an M-matrix.
 *
 *RETURN
 * -2 not well-constraint problem 
 * -1 inversion problem or no solution
 * 0  if not (or badly) midpoint-preconditioned 
 * 1  finish without problem
-------------------------------------------------- */

  int n = RowDimension (A);
  if (n != ColDimension (A))
    return -2;			// not well-constraint problem 

  MATRIX Id (n, n);
  MATRIX Zero (n, n);
  for (int i = 1; i <= n; i++)
    for (int j = 1; j <= n; j++)
      {
	Zero (i, j) = 0;
	Id (i, j) = i == j;
      }

  INTERVAL_MATRIX radA (A - Id);
  MATRIX InfA (Id - Abs (radA));
  MATRIX M (n, n);
  int res;
  if ((res = MatInverse (InfA, M)) != 1)
    return 0;

  if (!(M >= Zero))
    return 0;

  VECTOR b (Mid (B));
  VECTOR delta = Sup (B) - b;

  VECTOR xstar = M * (Abs (b) + delta);

  REAL xtildek, xutildek, nuk, max, min;

  for (int k = 1; k <= n; k++)
    {
      xtildek = (b (k) >= 0) ? xstar (k) : xstar (k) + 2 * M (k, k) * b (k);
      xutildek =
	(b (k) <= 0) ? -xstar (k) : -xstar (k) + 2 * M (k, k) * b (k);

      nuk = 1 / (2 * M (k, k) - 1);
      max = nuk * xutildek;
      if (max < 0)
	max = 0;

      min = nuk * xtildek;
      if (min > 0)
	min = 0;

      /* compute bounds of x(k) */
      if (xtildek >= max)
	{
	  if (xutildek <= min)
	    x (k) = INTERVAL (xutildek, xtildek);
	  else
	    x (k) = INTERVAL (max, xtildek);
	}
      else
	{
	  if (xutildek <= min)
	    x (k) = INTERVAL (xutildek, min);
	  else
	    {
	      return -1;
	    }
	}
    }
  return 1;
}

/*################################################
 * IV Tools
 *################################################*/

REAL
INorm (INTERVAL_VECTOR V)
{
  /* -------------------
   * Infinit norm of vector v
   ---------------------*/
  INT n = Dimension (V);
  REAL norme = Abs (V (1));
  REAL tmp;
  if (n > 1)
    for (int i = 2; i <= n; i++)
      {
	tmp = Abs (V (i));
	if (norme < tmp)
	  norme = tmp;
      }
  return norme;
}

REAL
INorm (INTERVAL_MATRIX M)
{
  /* -------------------
   * Infinit norm of matrix M
   ---------------------*/

  INT ln = RowDimension (M);
  INT cn = ColDimension (M);
  REAL norme, accu;
  // Particular cases
  if (ln == 1 & cn == 1)	// column number = row number = 1
    return Abs (M (1, 1));
  if (cn == 1)			// column number = 1
    return INorm (Col (M, 1));
  // General case
  norme = 0.0;
  for (int i = 1; i <= ln; i++)
    {
      accu = Abs (M (i, 1));
      for (int j = 2; j <= cn; j++)
	accu += Abs (M (i, j));
      if (accu > norme)
	norme = accu;
    }
  return norme;
}

int
MatInverse (MATRIX A, MATRIX & iA)
{
/* -------------------------
 * Inversion of a Matrix A in iA
 * Singularity return -1
 * else return 1
 --------------------------*/

  int code = 1;
  int wl, sl, old, flag_inverse;

  /* ----------------------------------
   * Errors management (Profile flags)
   -----------------------------------*/

  /* Stock previous values */
  wl = ErrorHandler::WarningLevel;
  sl = ErrorHandler::SevereLevel;
  old = ErrorHandler::LastErrorCode;

  /* Change value to avoid warnig */
  ErrorHandler::WarningLevel = 2500;
  ErrorHandler::SevereLevel = 3500;
  ErrorHandler::LastErrorCode = 0;

  /* Inverse Matrix */
  iA = Inverse (A);

  /* check the error */
  if (ErrorHandler::LastErrorCode == 1)
    code = -1;

  /* Put Error Handler codes in previous values */
  ErrorHandler::WarningLevel = wl;
  ErrorHandler::SevereLevel = sl;
  ErrorHandler::LastErrorCode = old;

  return code;

}

REAL
Perimeter (INTERVAL_VECTOR v)
{

  REAL p = Diam (v (1));

  INT k = Dimension (v);

  if (k > 1)
    for (int i = 2; i <= k; i++)
      p += Diam (v (i));

  return p;
}
