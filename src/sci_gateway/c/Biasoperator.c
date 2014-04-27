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

#include "Biasoperator.h"

#define Check(A,B,C); if ((A) == NULL) {free((B));free((C));*flag = -4; return -4;}
#define See(A,B); if ((A) == NULL) {free((B));*flag = -4; return -4;}

/*

// Flags convention for operations with Matrix intervals and Matrix reals :
// 2: operation concerning objects with different dimensions as matrix vector 
// multiplication or division
// 1: operation concerning objects with different dimensions including
// a one dimension object
// 0: operation concerning objects with similar dimensions
//-1: error flag, impossible operation
//-2: cowboy flag, there's something wrong !!
//-3: not implented yet ...
//-4: allocation problem

*/

/*

// Character conventions for binary operations :
// 'a' : addition
// 's' : substraction
// 'm' : multiplication
// 'x' : element by element multiplication
// 'k' : Kroneker multiplication
// 'r' : division
// 'd' : division element by element

// Character conventions for power fonctions :
// 'E' : element by element N-th power
// 'P' : N-th power
// 'R' :  element by element N-th root

// Character conventions for unary fonctions :
// '1' : absolute value
// '2' : exponential value
// '3' : neperian logarithmic value
// '4' : 10 th based logarithmic value
// 'a' : sinus
// 'b' : cosinus
// 'c' : tangente
// 'd' : cotangente
// 'e' : arcsinus
// 'f' : arccosinus
// 'g' : arctangente
// 'h' : arccotangente
// 'i' : hyperbolic sinus
// 'j' : hyperbolic cosinus
// 'k' : hyperbolic tangente
// 'l' : hyperbolic cotangente
// 'm' : hyperbolic arcsinus
// 'n' : hyperbolic arccosinus
// 'o' : hyperbolic arctangente
// 'p' : hyperbolic arccotangente

*/

int
power (const double *ainf, const double *asup, const int *sa,
       char op, const int *n, double *cinf, double *csup, int *flag)
{
  /*nth power of a interval matrix */
  if (op == 'P')
    {
      biaspowerNI4sci (ainf, asup, sa, n, cinf, csup, flag);
      return 1;
    }

  /* element by element n-th root */
  else if (op == 'R')
    {
      biasrootNI4sci (ainf, asup, sa, n, cinf, csup, flag);
      return 1;
    }

  /* element by element n-th power */
  else if (op == 'E')
    {
      biasRpowerNI4sci (ainf, asup, sa, n, cinf, csup, flag);
      return 1;
    }

  else
    {
      *flag = -2;
      return 1;
    }
}

int
unary (const double *ainf, const double *asup, const int *sa,
       char op, double *cinf, double *csup, int *flag)
{
  /* absolute value of a interval matrix */
  if (op == '1')
    {
      biasabsI4sci (ainf, asup, sa, cinf, csup, flag);
      return 1;
    }

  /* exponential value of a interval matrix */
  if (op == '2')
    {
      biasexpI4sci (ainf, asup, sa, cinf, csup, flag);
      return 1;
    }

  /* logarithmic value of a interval matrix */
  if (op == '3')
    {
      biaslogI4sci (ainf, asup, sa, cinf, csup, flag);
      return 1;
    }

  /* 10th based logarithmic value of a interval matrix */
  if (op == '4')
    {
      biaslog10I4sci (ainf, asup, sa, cinf, csup, flag);
      return 1;
    }

  /* sinus of a interval matrix */
  if (op == 'a')
    {
      biassinI4sci (ainf, asup, sa, cinf, csup, flag);
      return 1;
    }

  /* cosinus of a interval matrix */
  if (op == 'b')
    {
      biascosI4sci (ainf, asup, sa, cinf, csup, flag);
      return 1;
    }

  /* tangente of a interval matrix */
  if (op == 'c')
    {
      biastanI4sci (ainf, asup, sa, cinf, csup, flag);
      return 1;
    }

  /* cotangente of a interval matrix */
  if (op == 'd')
    {
      biascotI4sci (ainf, asup, sa, cinf, csup, flag);
      return 1;
    }

  /* arcsinus of a interval matrix */
  if (op == 'e')
    {
      biasarcsinI4sci (ainf, asup, sa, cinf, csup, flag);
      return 1;
    }

  /* arccosinus of a interval matrix */
  if (op == 'f')
    {
      biasarccosI4sci (ainf, asup, sa, cinf, csup, flag);
      return 1;
    }

  /* arctangente of a interval matrix */
  if (op == 'g')
    {
      biasarctanI4sci (ainf, asup, sa, cinf, csup, flag);
      return 1;
    }

  /* arccotangente of a interval matrix */
  if (op == 'h')
    {
      biasarccotI4sci (ainf, asup, sa, cinf, csup, flag);
      return 1;
    }

  /* hyperbolic sinus of a interval matrix */
  if (op == 'i')
    {
      biassinhI4sci (ainf, asup, sa, cinf, csup, flag);
      return 1;
    }

  /* hyperbolic cosinus of a interval matrix */
  if (op == 'j')
    {
      biascoshI4sci (ainf, asup, sa, cinf, csup, flag);
      return 1;
    }

  /* hyperbolic tangente of a interval matrix */
  if (op == 'k')
    {
      biastanhI4sci (ainf, asup, sa, cinf, csup, flag);
      return 1;
    }

  /* hyperbolic cotangente of a interval matrix */
  if (op == 'l')
    {
      biascothI4sci (ainf, asup, sa, cinf, csup, flag);
      return 1;
    }

  /* hyperbolic arcsinus of a interval matrix */
  if (op == 'm')
    {
      biasarsinhI4sci (ainf, asup, sa, cinf, csup, flag);
      return 1;
    }

  /* hyperbolic arccosinus of a interval matrix */
  if (op == 'n')
    {
      biasarcoshI4sci (ainf, asup, sa, cinf, csup, flag);
      return 1;
    }

  /* hyperbolic arctangente of a interval matrix */
  if (op == 'o')
    {
      biasartanhI4sci (ainf, asup, sa, cinf, csup, flag);
      return 1;
    }

  /* hyperbolic arccotangente of a interval matrix */
  if (op == 'p')
    {
      biasarcothI4sci (ainf, asup, sa, cinf, csup, flag);
      return 1;
    }

}

int
operate (const double *ainf, const double *asup, const int *sa,
	 const double *binf, const double *bsup, const int *sb,
	 char op, double *cinf, double *csup, int *flag)
{
  if (op == 'a')
    {
      /* AddII */
      if (bsup != NULL)
	{
	  biasaddII4sci (ainf, asup, sa, binf, bsup, sb, cinf, csup, flag);
	  return 1;
	}

      /* AddIR & AddRI */
      else
	{
	  biasaddIR4sci (ainf, asup, sa, binf, sb, cinf, csup, flag);
	  return 1;
	}
    }
  else if (op == 's')
    {
      /* SubII */
      if (bsup != NULL)
	{
	  biassubII4sci (ainf, asup, sa, binf, bsup, sb, cinf, csup, flag);
	  return 1;
	}

      /* SubIR & SubRI */
      else
	{
	  biassubIR4sci (ainf, asup, sa, binf, sb, cinf, csup, flag);
	  return 1;
	}
    }
  else if (op == 'm')
    {
      /* MulII */
      if (bsup != NULL && asup != NULL)
	{
	  biasmulII4sci (ainf, asup, sa, binf, bsup, sb, cinf, csup, flag);
	  return 1;
	}

      /* MulIR */
      else if (bsup == NULL)
	{
	  biasmulIR4sci (ainf, asup, sa, binf, sb, cinf, csup, flag);
	  return 1;
	}

      /* MulRI */
      else if (asup == NULL)
	{
	  biasmulRI4sci (ainf, sa, binf, bsup, sb, cinf, csup, flag);
	  return 1;
	}

      else
	{
	  *flag = -2;
	  return 1;
	}
    }

  else if (op == 'x')
    {
      /* MulII element by element */
      if (bsup != NULL)
	{
	  biasRmulII4sci (ainf, asup, sa, binf, bsup, sb, cinf, csup, flag);
	  return 1;
	}

      /* MulIR element by element */
      else
	{
	  biasRmulIR4sci (ainf, asup, sa, binf, sb, cinf, csup, flag);
	  return 1;
	}
    }

  else if (op == 'k')
    {
      /* Kroneker MulII */
      if (bsup != NULL && asup != NULL)
	{
	  biasRLmulII4sci (ainf, asup, sa, binf, bsup, sb, cinf, csup, flag);
	  return 1;
	}

      /* Kroneker MulIR */
      else if (bsup == NULL)
	{
	  biasRLmulIR4sci (ainf, asup, sa, binf, sb, cinf, csup, flag);
	  return 1;
	}

      /* Kroneker MulRI */
      else if (asup == NULL)
	{
	  biasRLmulRI4sci (ainf, sa, binf, bsup, sb, cinf, csup, flag);
	  return 1;
	}

      else
	{
	  *flag = -2;
	  return 1;
	}
    }

  else if (op == 'r')
    {
      /* Division of two intervals */
      if (bsup != NULL && asup != NULL)
	{
	  biasdivII4sci (ainf, asup, sa, binf, bsup, sb, cinf, csup, flag);
	  return 1;
	}

      /* Division of an interval by a real */
      else if (bsup == NULL)
	{
	  biasdivIR4sci (ainf, asup, sa, binf, sb, cinf, csup, flag);
	  return 1;
	}

      /*Division of a real by an interval */
      else if (asup == NULL)
	{
	  biasdivRI4sci (ainf, sa, binf, bsup, sb, cinf, csup, flag);
	  return 1;
	}
      else
	{
	  *flag = -2;
	  return 1;
	}
    }

  else if (op == 'd')
    {
      /* Division element by element of two intervals */
      if (bsup != NULL && asup != NULL)
	{
	  biasRdivII4sci (ainf, asup, sa, binf, bsup, sb, cinf, csup, flag);
	  return 1;
	}

      /* Division element by element of an interval by a real */
      else if (bsup == NULL)
	{
	  biasRdivIR4sci (ainf, asup, sa, binf, sb, cinf, csup, flag);
	  return 1;
	}

      /*Division element by element of a real by an interval */
      else if (asup == NULL)
	{
	  biasRdivRI4sci (ainf, sa, binf, bsup, sb, cinf, csup, flag);
	  return 1;
	}
      else
	{
	  *flag = -2;
	  return 1;
	}
    }
  else
    {
      *flag = -2;
      return 1;
    }
}

int
biasaddII4sci (const double *ainf, const double *asup, const int *sa,
	       const double *binf, const double *bsup, const int *sb,
	       double *cinf, double *csup, int *flag)
{
  //feclearexcept (FE_ALL_EXCEPT);
  int i;
  BIASINTERVAL *pA;
  BIASINTERVAL *pB;
  BIASINTERVAL *pC;

  /* Case with different dimensions */
  if (sa[0] != sb[0] || sa[1] != sb[1])
    {
      /*add a single interval */

      if (sb[0] == 1 && sb[1] == 1)
	{
	  pA =
	    (BIASINTERVAL *) malloc (sa[0] * sa[1] * sizeof (BIASINTERVAL));
	  pB =
	    (BIASINTERVAL *) malloc (sb[0] * sb[1] * sizeof (BIASINTERVAL));
	  pC =
	    (BIASINTERVAL *) malloc (max (sa[0] * sa[1], sb[0] * sb[1]) *
				     sizeof (BIASINTERVAL));
	  Check (pA, pB, pC);
	  Check (pB, pA, pC);
	  Check (pC, pA, pB);

	  pB[0].inf = binf[0];
	  pB[0].sup = bsup[0];

	  for (i = 0; i < sa[0] * sa[1]; i++)
	    {
	      pA[i].inf = ainf[i];
	      pA[i].sup = asup[i];
	      BiasAddII (&(pC[i]), &(pA[i]), &(pB[0]));
	      cinf[i] = pC[i].inf;
	      csup[i] = pC[i].sup;
	    }

	  free (pA);
	  free (pB);
	  free (pC);
	  *flag = 1;
	  return 1;
	}
      else if (sa[0] == 1 && sa[1] == 1)
	{
	  pA =
	    (BIASINTERVAL *) malloc (sa[0] * sa[1] * sizeof (BIASINTERVAL));
	  pB =
	    (BIASINTERVAL *) malloc (sb[0] * sb[1] * sizeof (BIASINTERVAL));
	  pC =
	    (BIASINTERVAL *) malloc (max (sa[0] * sa[1], sb[0] * sb[1]) *
				     sizeof (BIASINTERVAL));
	  Check (pA, pB, pC);
	  Check (pB, pA, pC);
	  Check (pC, pA, pB);

	  pA[0].inf = ainf[0];
	  pA[0].sup = asup[0];

	  for (i = 0; i < sb[0] * sb[1]; i++)
	    {
	      pB[i].inf = binf[i];
	      pB[i].sup = bsup[i];
	      BiasAddII (&(pC[i]), &(pA[0]), &(pB[i]));
	      cinf[i] = pC[i].inf;
	      csup[i] = pC[i].sup;
	    }

	  free (pA);
	  free (pB);
	  free (pC);
	  *flag = 1;
	  return 1;
	}

      /* Impossible operation */

      else
	{
	  //printf (" INTERVAL Dimension error between matrix");
	  *flag = -1;
	  return -1;
	}
    }

  /*Case with similar dimensions */
  else
    {
      pA = (BIASINTERVAL *) malloc (sa[0] * sa[1] * sizeof (BIASINTERVAL));
      pB = (BIASINTERVAL *) malloc (sb[0] * sb[1] * sizeof (BIASINTERVAL));
      pC =
	(BIASINTERVAL *) malloc (max (sa[0] * sa[1], sb[0] * sb[1]) *
				 sizeof (BIASINTERVAL));
      Check (pA, pB, pC);
      Check (pB, pA, pC);
      Check (pC, pA, pB);

      for (i = 0; i < sa[0] * sa[1]; i++)
	{
	  pA[i].inf = ainf[i];
	  pA[i].sup = asup[i];
	  pB[i].inf = binf[i];
	  pB[i].sup = bsup[i];
	  pC[i].inf = 1;
	  pC[i].sup = 1;
	}

      BiasAddMIMI (pC, pA, pB, sa[0], sa[1]);

      for (i = 0; i < sa[0] * sa[1]; i++)
	{
	  cinf[i] = pC[i].inf;
	  csup[i] = pC[i].sup;
	}

      free (pA);
      free (pB);
      free (pC);
      *flag = 0;
      return 0;
    }
}



int
biasaddIR4sci (const double *ainf, const double *asup, const int *sa,
	       const double *b, const int *sb,
	       double *cinf, double *csup, int *flag)
{
  //feclearexcept (FE_ALL_EXCEPT);
  int i;
  BIASINTERVAL *pA;
  BIASINTERVAL *pC;

  /* Case with different dimensions */
  if (sa[0] != sb[0] || sa[1] != sb[1])
    {
      /*add a single real */

      if (sb[0] == 1 && sb[1] == 1)
	{
	  pA =
	    (BIASINTERVAL *) malloc (sa[0] * sa[1] * sizeof (BIASINTERVAL));
	  pC =
	    (BIASINTERVAL *) malloc (max (sa[0] * sa[1], sb[0] * sb[1]) *
				     sizeof (BIASINTERVAL));
	  See (pA, pC);
	  See (pC, pA);

	  for (i = 0; i < sa[0] * sa[1]; i++)
	    {
	      pA[i].inf = ainf[i];
	      pA[i].sup = asup[i];
	      BiasAddIR (&(pC[i]), &(pA[i]), &(b[0]));
	      cinf[i] = pC[i].inf;
	      csup[i] = pC[i].sup;
	    }

	  free (pA);
	  free (pC);
	  *flag = 1;
	  return 1;
	}
      /*add a single interval */

      else if (sa[0] == 1 && sa[1] == 1)
	{
	  pA =
	    (BIASINTERVAL *) malloc (sa[0] * sa[1] * sizeof (BIASINTERVAL));
	  pC =
	    (BIASINTERVAL *) malloc (max (sa[0] * sa[1], sb[0] * sb[1]) *
				     sizeof (BIASINTERVAL));
	  See (pA, pC);
	  See (pC, pA);

	  pA[0].inf = ainf[0];
	  pA[0].sup = asup[0];

	  for (i = 0; i < sb[0] * sb[1]; i++)
	    {
	      BiasAddIR (&(pC[i]), &(pA[0]), &(b[i]));
	      cinf[i] = pC[i].inf;
	      csup[i] = pC[i].sup;
	    }

	  free (pA);
	  free (pC);
	  *flag = 1;
	  return 1;
	}

      /*Impossible operation */
      else
	{
	  //printf (" INTERVAL Dimension error between matrix");
	  *flag = -1;
	  return -1;
	}
    }

  /*Case with similar dimensions */
  else
    {
      pA = (BIASINTERVAL *) malloc (sa[0] * sa[1] * sizeof (BIASINTERVAL));
      pC =
	(BIASINTERVAL *) malloc (max (sa[0] * sa[1], sb[0] * sb[1]) *
				 sizeof (BIASINTERVAL));
      See (pA, pC);
      See (pC, pA);

      for (i = 0; i < sa[0] * sa[1]; i++)
	{
	  pA[i].inf = ainf[i];
	  pA[i].sup = asup[i];
	}

      BiasAddMIMR (pC, pA, b, sa[0], sa[1]);

      for (i = 0; i < sa[0] * sa[1]; i++)
	{
	  cinf[i] = pC[i].inf;
	  csup[i] = pC[i].sup;
	}

      free (pA);
      free (pC);
      *flag = 0;
      return 0;
    }

}



int
biassubII4sci (const double *ainf, const double *asup, const int *sa,
	       const double *binf, const double *bsup, const int *sb,
	       double *cinf, double *csup, int *flag)
{
  //feclearexcept (FE_ALL_EXCEPT);
  int i;
  BIASINTERVAL *pA;
  BIASINTERVAL *pB;
  BIASINTERVAL *pC;

  /* Case with different dimensions */
  if (sa[0] != sb[0] || sa[1] != sb[1])
    {
      /*sub a single interval */

      if (sb[0] == 1 && sb[1] == 1)
	{
	  pA =
	    (BIASINTERVAL *) malloc (sa[0] * sa[1] * sizeof (BIASINTERVAL));
	  pB =
	    (BIASINTERVAL *) malloc (sb[0] * sb[1] * sizeof (BIASINTERVAL));
	  pC =
	    (BIASINTERVAL *) malloc (max (sa[0] * sa[1], sb[0] * sb[1]) *
				     sizeof (BIASINTERVAL));
	  Check (pA, pB, pC);
	  Check (pB, pA, pC);
	  Check (pC, pA, pB);

	  pB[0].inf = binf[0];
	  pB[0].sup = bsup[0];

	  for (i = 0; i < sa[0] * sa[1]; i++)
	    {
	      pA[i].inf = ainf[i];
	      pA[i].sup = asup[i];
	      BiasSubII (&(pC[i]), &(pA[i]), &(pB[0]));
	      cinf[i] = pC[i].inf;
	      csup[i] = pC[i].sup;
	    }

	  free (pA);
	  free (pB);
	  free (pC);
	  *flag = 1;
	  return 1;
	}
      else if (sa[0] == 1 && sa[1] == 1)
	{
	  pA =
	    (BIASINTERVAL *) malloc (sa[0] * sa[1] * sizeof (BIASINTERVAL));
	  pB =
	    (BIASINTERVAL *) malloc (sb[0] * sb[1] * sizeof (BIASINTERVAL));
	  pC =
	    (BIASINTERVAL *) malloc (max (sa[0] * sa[1], sb[0] * sb[1]) *
				     sizeof (BIASINTERVAL));
	  Check (pA, pB, pC);
	  Check (pB, pA, pC);
	  Check (pC, pA, pB);

	  pA[0].inf = ainf[0];
	  pA[0].sup = asup[0];

	  for (i = 0; i < sb[0] * sb[1]; i++)
	    {
	      pB[i].inf = binf[i];
	      pB[i].sup = bsup[i];
	      BiasSubII (&(pC[i]), &(pA[0]), &(pB[i]));
	      cinf[i] = pC[i].inf;
	      csup[i] = pC[i].sup;
	    }

	  free (pA);
	  free (pB);
	  free (pC);
	  *flag = 1;
	  return 1;
	}

      /* Impossible operation */

      else
	{
	  //printf (" INTERVAL Dimension error between matrix");
	  *flag = -1;
	  return -1;
	}
    }

  /*Case with similar dimensions */
  else
    {
      pA = (BIASINTERVAL *) malloc (sa[0] * sa[1] * sizeof (BIASINTERVAL));
      pB = (BIASINTERVAL *) malloc (sb[0] * sb[1] * sizeof (BIASINTERVAL));
      pC =
	(BIASINTERVAL *) malloc (max (sa[0] * sa[1], sb[0] * sb[1]) *
				 sizeof (BIASINTERVAL));
      Check (pA, pB, pC);
      Check (pB, pA, pC);
      Check (pC, pA, pB);

      for (i = 0; i < sa[0] * sa[1]; i++)
	{
	  pA[i].inf = ainf[i];
	  pA[i].sup = asup[i];
	  pB[i].inf = binf[i];
	  pB[i].sup = bsup[i];
	}

      BiasSubMIMI (pC, pA, pB, sa[0], sa[1]);

      for (i = 0; i < sa[0] * sa[1]; i++)
	{
	  cinf[i] = pC[i].inf;
	  csup[i] = pC[i].sup;
	}

      free (pA);
      free (pB);
      free (pC);
      *flag = 0;
      return 0;
    }
}



int
biassubIR4sci (const double *ainf, const double *asup, const int *sa,
	       const double *b, const int *sb,
	       double *cinf, double *csup, int *flag)
{
  //feclearexcept (FE_ALL_EXCEPT);
  int i;
  BIASINTERVAL *pA;
  BIASINTERVAL *pC;

  /* Case with different dimensions */
  if (sa[0] != sb[0] || sa[1] != sb[1])
    {
      /*sub a single real */

      if (sb[0] == 1 && sb[1] == 1)
	{
	  pA =
	    (BIASINTERVAL *) malloc (sa[0] * sa[1] * sizeof (BIASINTERVAL));
	  pC =
	    (BIASINTERVAL *) malloc (max (sa[0] * sa[1], sb[0] * sb[1]) *
				     sizeof (BIASINTERVAL));
	  See (pA, pC);
	  See (pC, pA);

	  for (i = 0; i < sa[0] * sa[1]; i++)
	    {
	      pA[i].inf = ainf[i];
	      pA[i].sup = asup[i];
	      BiasSubIR (&(pC[i]), &(pA[i]), &(b[0]));
	      cinf[i] = pC[i].inf;
	      csup[i] = pC[i].sup;
	    }

	  free (pA);
	  free (pC);
	  *flag = 1;
	  return 1;
	}
      /*sub a single interval */

      else if (sa[0] == 1 && sa[1] == 1)
	{
	  pA =
	    (BIASINTERVAL *) malloc (sa[0] * sa[1] * sizeof (BIASINTERVAL));
	  pC =
	    (BIASINTERVAL *) malloc (max (sa[0] * sa[1], sb[0] * sb[1]) *
				     sizeof (BIASINTERVAL));
	  See (pA, pC);
	  See (pC, pA);

	  pA[0].inf = ainf[0];
	  pA[0].sup = asup[0];

	  for (i = 0; i < sb[0] * sb[1]; i++)
	    {
	      BiasSubIR (&(pC[i]), &(pA[0]), &(b[i]));
	      cinf[i] = pC[i].inf;
	      csup[i] = pC[i].sup;
	    }

	  free (pA);
	  free (pC);
	  *flag = 1;
	  return 1;
	}

      /*Impossible operation */
      else
	{
	  //printf (" INTERVAL Dimension error between matrix");
	  *flag = -1;
	  return -1;
	}
    }

  /*Case with similar dimensions */
  else
    {
      pA = (BIASINTERVAL *) malloc (sa[0] * sa[1] * sizeof (BIASINTERVAL));
      pC =
	(BIASINTERVAL *) malloc (max (sa[0] * sa[1], sb[0] * sb[1]) *
				 sizeof (BIASINTERVAL));
      See (pA, pC);
      See (pC, pA);

      for (i = 0; i < sa[0] * sa[1]; i++)
	{
	  pA[i].inf = ainf[i];
	  pA[i].sup = asup[i];
	}

      BiasSubMIMR (pC, pA, b, sa[0], sa[1]);

      for (i = 0; i < sa[0] * sa[1]; i++)
	{
	  cinf[i] = pC[i].inf;
	  csup[i] = pC[i].sup;
	}

      free (pA);
      free (pC);
      *flag = 0;
      return 0;
    }
}


int
biasmulII4sci (const double *ainf, const double *asup, const int *sa,
	       const double *binf, const double *bsup, const int *sb,
	       double *cinf, double *csup, int *flag)
{
  //feclearexcept (FE_ALL_EXCEPT);

  int i;
  BIASINTERVAL *pA;
  BIASINTERVAL *pB;
  BIASINTERVAL *pC;

  /*Invalid Bias operation control constants */
  int rinf = 0;
  int rsup = 0;

  /* Case with dimensions not conform for
     matrix multiplication */

  if (sa[1] != sb[0])
    {
      /*mul a single interval */

      if (sb[0] == 1 && sb[1] == 1)
	{
	  pA =
	    (BIASINTERVAL *) malloc (sa[0] * sa[1] * sizeof (BIASINTERVAL));
	  pB =
	    (BIASINTERVAL *) malloc (sb[0] * sb[1] * sizeof (BIASINTERVAL));
	  pC =
	    (BIASINTERVAL *) malloc (max (sa[0] * sa[1], sb[0] * sb[1]) *
				     sizeof (BIASINTERVAL));
	  Check (pA, pB, pC);
	  Check (pB, pA, pC);
	  Check (pC, pA, pB);

	  pB[0].inf = binf[0];
	  pB[0].sup = bsup[0];

	  /*Invalid Bias operation control */

	  if (pB[0].inf == 0 || fabs (pB[0].inf) == HUGE_VAL)
	    {
	      rinf++;
	    }

	  if (pB[0].sup == 0 || fabs (pB[0].sup) == HUGE_VAL)
	    {
	      rsup++;
	    }

	  for (i = 0; i < sa[0] * sa[1]; i++)
	    {
	      pA[i].inf = ainf[iterf2c (i, sa[0], sa[1])];
	      pA[i].sup = asup[iterf2c (i, sa[0], sa[1])];

	      /*Invalid Bias operation control */

	      if (rinf != 0 || rsup != 0)
		{
		  if (pA[i].inf == 0 || fabs (pA[i].inf) == HUGE_VAL)
		    {
		      rinf++;
		    }

		  if (pA[i].sup == 0 || fabs (pA[i].sup) == HUGE_VAL)
		    {
		      rsup++;
		    }
		}
	    }

	  if (rinf >= 1 && rsup >= 1)
	    {
	      /*Case with invalid Bias operations */

	      DIYMulIMI (pC, &(pB[0]), pA, sa);

	      for (i = 0; i < sa[0] * sa[1]; i++)
		{
		  cinf[iterf2c (i, sa[0], sa[1])] = pC[i].inf;
		  csup[iterf2c (i, sa[0], sa[1])] = pC[i].sup;

		}

	      free (pA);
	      free (pB);
	      free (pC);
	      *flag = 1;
	      return 1;
	    }
	  else
	    {
	      /*General valid case */

	      BiasMulIMI (pC, &(pB[0]), pA, sa[0], sa[1]);

	      for (i = 0; i < sa[0] * sa[1]; i++)
		{
		  cinf[iterf2c (i, sa[0], sa[1])] = pC[i].inf;
		  csup[iterf2c (i, sa[0], sa[1])] = pC[i].sup;

		}

	      free (pA);
	      free (pB);
	      free (pC);
	      *flag = 1;
	      return 1;
	    }
	}

      else if (sa[0] == 1 && sa[1] == 1)
	{
	  pA =
	    (BIASINTERVAL *) malloc (sa[0] * sa[1] * sizeof (BIASINTERVAL));
	  pB =
	    (BIASINTERVAL *) malloc (sb[0] * sb[1] * sizeof (BIASINTERVAL));
	  pC =
	    (BIASINTERVAL *) malloc (max (sa[0] * sa[1], sb[0] * sb[1]) *
				     sizeof (BIASINTERVAL));
	  Check (pA, pB, pC);
	  Check (pB, pA, pC);
	  Check (pC, pA, pB);

	  pA[0].inf = ainf[0];
	  pA[0].sup = asup[0];

	  /*Invalid Bias operation control */

	  if (pA[0].inf == 0 || fabs (pA[0].inf) == HUGE_VAL)
	    {
	      rinf++;
	    }

	  if (pA[0].sup == 0 || fabs (pA[0].sup) == HUGE_VAL)
	    {
	      rsup++;
	    }

	  for (i = 0; i < sb[0] * sb[1]; i++)
	    {
	      pB[i].inf = binf[iterf2c (i, sb[0], sb[1])];
	      pB[i].sup = bsup[iterf2c (i, sb[0], sb[1])];

	      /*Invalid Bias operation control */

	      if (rinf != 0 || rsup != 0)
		{
		  if (pB[i].inf == 0 || fabs (pB[i].inf) == HUGE_VAL)
		    {
		      rinf++;
		    }

		  if (pB[i].sup == 0 || fabs (pB[i].sup) == HUGE_VAL)
		    {
		      rsup++;
		    }
		}
	    }

	  if (rinf >= 1 && rsup >= 1)
	    {
	      /*Case with invalid Bias operations */

	      DIYMulIMI (pC, &(pA[0]), pB, sa);

	      for (i = 0; i < sa[0] * sa[1]; i++)
		{
		  cinf[iterf2c (i, sa[0], sa[1])] = pC[i].inf;
		  csup[iterf2c (i, sa[0], sa[1])] = pC[i].sup;
		}

	      free (pA);
	      free (pB);
	      free (pC);
	      *flag = 1;
	      return 1;
	    }
	  else
	    {
	      /*General valid case */

	      BiasMulIMI (pC, &(pA[0]), pB, sb[0], sb[1]);

	      for (i = 0; i < sb[0] * sb[1]; i++)
		{
		  cinf[iterf2c (i, sb[0], sb[1])] = pC[i].inf;
		  csup[iterf2c (i, sb[0], sb[1])] = pC[i].sup;
		}

	      free (pA);
	      free (pB);
	      free (pC);
	      *flag = 1;
	      return 1;
	    }
	}

      /* Impossible operation */

      else
	{
	  //printf (" INTERVAL Dimension error between matrix");
	  *flag = -1;
	  return -1;
	}
    }

  /*Case with conform dimensions */
  else
    {
      pA = (BIASINTERVAL *) malloc (sa[0] * sa[1] * sizeof (BIASINTERVAL));
      pB = (BIASINTERVAL *) malloc (sb[0] * sb[1] * sizeof (BIASINTERVAL));
      pC = (BIASINTERVAL *) malloc (sa[0] * sb[1] * sizeof (BIASINTERVAL));
      Check (pA, pB, pC);
      Check (pB, pA, pC);
      Check (pC, pA, pB);

      for (i = 0; i < sa[0] * sa[1]; i++)
	{
	  pA[i].inf = ainf[iterf2c (i, sa[0], sa[1])];
	  pA[i].sup = asup[iterf2c (i, sa[0], sa[1])];

	  /*Invalid Bias operation control */

	  if (pA[i].inf == 0 || fabs (pA[i].inf) == HUGE_VAL)
	    {
	      rinf++;
	    }

	  if (pA[i].sup == 0 || fabs (pA[i].sup) == HUGE_VAL)
	    {
	      rsup++;
	    }
	}

      for (i = 0; i < sb[0] * sb[1]; i++)
	{
	  pB[i].inf = binf[iterf2c (i, sb[0], sb[1])];
	  pB[i].sup = bsup[iterf2c (i, sb[0], sb[1])];

	  /*Invalid Bias operation control */

	  if (rinf != 0 || rsup != 0)
	    {
	      if (pB[i].inf == 0 || fabs (pB[i].inf) == HUGE_VAL)
		{
		  rinf++;
		}

	      if (pB[i].sup == 0 || fabs (pB[i].sup) == HUGE_VAL)
		{
		  rsup++;
		}
	    }
	}

      if (rinf >= 1 && rsup >= 1)
	{
	  /*Case with invalid Bias operations */

	  DIYMulMIMI (pC, pA, pB, sa, sb);

	  for (i = 0; i < sa[0] * sb[1]; i++)
	    {
	      cinf[iterf2c (i, sa[0], sb[1])] = pC[i].inf;
	      csup[iterf2c (i, sa[0], sb[1])] = pC[i].sup;
	    }

	  free (pA);
	  free (pB);
	  free (pC);
	  *flag = 2;
	  return 0;
	}
      else
	{
	  /*General valid case */

	  BiasMulMIMI (pC, pA, pB, sa[0], sa[1], sb[1]);

	  for (i = 0; i < sa[0] * sb[1]; i++)
	    {
	      cinf[iterf2c (i, sa[0], sb[1])] = pC[i].inf;
	      csup[iterf2c (i, sa[0], sb[1])] = pC[i].sup;
	    }

	  free (pA);
	  free (pB);
	  free (pC);
	  *flag = 2;
	  return 0;
	}
    }
}

int
biasmulIR4sci (const double *ainf, const double *asup, const int *sa,
	       const double *b, const int *sb,
	       double *cinf, double *csup, int *flag)
{
  //feclearexcept (FE_ALL_EXCEPT);
  int i;
  BIASINTERVAL *pA;
  BIASINTERVAL *pC;
  double *rB;

  /*Invalid Bias operation control constants */
  int rinf = 0;
  int rsup = 0;

  /* Case with dimensions not conform for
     matrix multiplication */

  if (sa[1] != sb[0])
    {
      /*mul a single real */

      if (sb[0] == 1 && sb[1] == 1)
	{
	  /*Invalid Bias operation control */

	  if (b[0] == 0 || fabs (b[0]) == HUGE_VAL)
	    {
	      rinf++;
	      rsup++;
	    }

	  pA =
	    (BIASINTERVAL *) malloc (sa[0] * sa[1] * sizeof (BIASINTERVAL));
	  pC =
	    (BIASINTERVAL *) malloc (max (sa[0] * sa[1], sb[0] * sb[1]) *
				     sizeof (BIASINTERVAL));
	  See (pA, pC);
	  See (pC, pA);

	  for (i = 0; i < sa[0] * sa[1]; i++)
	    {
	      pA[i].inf = ainf[iterf2c (i, sa[0], sa[1])];
	      pA[i].sup = asup[iterf2c (i, sa[0], sa[1])];

	      /*Invalid Bias operation control */

	      if (rinf != 0 || rsup != 0)
		{
		  if (pA[i].inf == 0 || fabs (pA[i].inf) == HUGE_VAL)
		    {
		      rinf++;
		    }

		  if (pA[i].sup == 0 || fabs (pA[i].sup) == HUGE_VAL)
		    {
		      rsup++;
		    }
		}
	    }

	  if (rinf > 1 && rsup > 1)
	    {
	      /*Case with invalid Bias operations */

	      DIYMulRMI (pC, &(b[0]), pA, sa);

	      for (i = 0; i < sa[0] * sa[1]; i++)
		{
		  cinf[iterf2c (i, sa[0], sa[1])] = pC[i].inf;
		  csup[iterf2c (i, sa[0], sa[1])] = pC[i].sup;
		}

	      free (pA);
	      free (pC);
	      *flag = 1;
	      return 1;
	    }
	  else
	    {
	      /*General valid case */

	      BiasMulRMI (pC, &(b[0]), pA, sa[0], sa[1]);

	      for (i = 0; i < sa[0] * sa[1]; i++)
		{
		  cinf[iterf2c (i, sa[0], sa[1])] = pC[i].inf;
		  csup[iterf2c (i, sa[0], sa[1])] = pC[i].sup;
		}

	      free (pA);
	      free (pC);
	      *flag = 1;
	      return 1;
	    }
	}
      /*mul a single interval */

      else if (sa[0] == 1 && sa[1] == 1)
	{
	  pA =
	    (BIASINTERVAL *) malloc (sa[0] * sa[1] * sizeof (BIASINTERVAL));
	  pC =
	    (BIASINTERVAL *) malloc (max (sa[0] * sa[1], sb[0] * sb[1]) *
				     sizeof (BIASINTERVAL));
	  rB = (double *) malloc (sb[0] * sb[1] * sizeof (double));
	  Check (pA, rB, pC);
	  Check (rB, pA, pC);
	  Check (pC, pA, rB);

	  pA[0].inf = ainf[0];
	  pA[0].sup = asup[0];

	  /*Invalid Bias operation control */

	  if (pA[0].inf == 0 || fabs (pA[0].inf) == HUGE_VAL)
	    {
	      rinf++;
	    }

	  if (pA[0].sup == 0 || fabs (pA[0].sup) == HUGE_VAL)
	    {
	      rsup++;
	    }

	  for (i = 0; i < sb[0] * sb[1]; i++)
	    {
	      if (rinf != 0 || rsup != 0)
		{
		  if (b[i] == 0 || fabs (b[i]) == HUGE_VAL)
		    {
		      rinf++;
		      rsup++;
		    }
		}
	    }

	  if (rinf >= 2 && rsup >= 2)
	    {
	      /*Case with invalid Bias operations */

	      DIYMulIMR (pC, &(pA[0]), b, sb);

	      for (i = 0; i < sb[0] * sb[1]; i++)
		{
		  cinf[iterf2c (i, sb[0], sb[1])] = pC[i].inf;
		  csup[iterf2c (i, sb[0], sb[1])] = pC[i].sup;
		}

	      free (pA);
	      free (pC);
	      free (rB);
	      *flag = 1;
	      return 1;
	    }

	  else
	    {
	      /*General valid case */

	      BiasMulIMR (pC, &(pA[0]), b, sb[0], sb[1]);

	      for (i = 0; i < sb[0] * sb[1]; i++)
		{
		  cinf[iterf2c (i, sb[0], sb[1])] = pC[i].inf;
		  csup[iterf2c (i, sb[0], sb[1])] = pC[i].sup;
		}

	      free (pA);
	      free (pC);
	      free (rB);
	      *flag = 1;
	      return 1;
	    }
	}

      /*Impossible operation */
      else
	{
	  //printf (" INTERVAL Dimension error between matrix");
	  *flag = -1;
	  return -1;
	}
    }

  /*Case with conform dimensions */
  else
    {
      pA = (BIASINTERVAL *) malloc (sa[0] * sa[1] * sizeof (BIASINTERVAL));
      pC = (BIASINTERVAL *) malloc (sa[0] * sb[1] * sizeof (BIASINTERVAL));
      rB = (double *) malloc (sb[0] * sb[1] * sizeof (double));
      Check (pA, rB, pC);
      Check (rB, pA, pC);
      Check (pC, pA, rB);

      for (i = 0; i < sa[0] * sa[1]; i++)
	{
	  pA[i].inf = ainf[iterf2c (i, sa[0], sa[1])];
	  pA[i].sup = asup[iterf2c (i, sa[0], sa[1])];

	  /*Invalid Bias operation control */

	  if (pA[i].inf == 0 || fabs (pA[i].inf) == HUGE_VAL)
	    {
	      rinf++;
	    }

	  if (pA[i].sup == 0 || fabs (pA[i].sup) == HUGE_VAL)
	    {
	      rsup++;
	    }
	}

      for (i = 0; i < sb[0] * sb[1]; i++)
	{
	  rB[i] = b[iterf2c (i, sb[0], sb[1])];

	  /*Invalid Bias operation control */

	  if (rinf != 0 || rsup != 0)
	    {
	      if (rB[i] == 0 || fabs (rB[i]) == HUGE_VAL)
		{
		  rinf++;
		  rsup++;
		}
	    }
	}

      if (rinf >= 1 && rsup >= 1)
	{
	  /*Case with invalid Bias operations */

	  DIYMulMIMR (pC, pA, rB, sa, sb);

	  for (i = 0; i < sa[0] * sb[1]; i++)
	    {
	      cinf[iterf2c (i, sa[0], sb[1])] = pC[i].inf;
	      csup[iterf2c (i, sa[0], sb[1])] = pC[i].sup;
	    }

	  free (pA);
	  free (pC);
	  free (rB);
	  *flag = 2;
	  return 0;
	}

      else
	{
	  /*General valid case */

	  BiasMulMIMR (pC, pA, rB, sa[0], sa[1], sb[1]);

	  for (i = 0; i < sa[0] * sb[1]; i++)
	    {
	      cinf[iterf2c (i, sa[0], sb[1])] = pC[i].inf;
	      csup[iterf2c (i, sa[0], sb[1])] = pC[i].sup;
	    }

	  free (pA);
	  free (pC);
	  free (rB);
	  *flag = 2;
	  return 0;
	}
    }
}


int
biasmulRI4sci (const double *a, const int *sa,
	       const double *binf, const double *bsup, const int *sb,
	       double *cinf, double *csup, int *flag)
{
  //feclearexcept (FE_ALL_EXCEPT);
  int i;
  BIASINTERVAL *pB;
  BIASINTERVAL *pC;
  double *rA;

  /*Invalid Bias operation control constants */
  int rinf = 0;
  int rsup = 0;

  /* Case with dimensions not conform for
     matrix multiplication */

  if (sa[1] != sb[0])
    {

      biasmulIR4sci (binf, bsup, sb, a, sa, cinf, csup, flag);
      return 0;

    }

  /*Case with conform dimensions */
  else
    {
      pB = (BIASINTERVAL *) malloc (sb[0] * sb[1] * sizeof (BIASINTERVAL));
      pC = (BIASINTERVAL *) malloc (sa[0] * sb[1] * sizeof (BIASINTERVAL));
      rA = (double *) malloc (sa[0] * sa[1] * sizeof (double));
      Check (rA, pB, pC);
      Check (pB, rA, pC);
      Check (pC, rA, pB);

      for (i = 0; i < sb[0] * sb[1]; i++)
	{
	  pB[i].inf = binf[iterf2c (i, sb[0], sb[1])];
	  pB[i].sup = bsup[iterf2c (i, sb[0], sb[1])];

	  /*Invalid Bias operation control */

	  if (pB[i].inf == 0 || fabs (pB[i].inf) == HUGE_VAL)
	    {
	      rinf++;
	    }

	  if (pB[i].sup == 0 || fabs (pB[i].sup) == HUGE_VAL)
	    {
	      rsup++;
	    }
	}

      for (i = 0; i < sa[0] * sa[1]; i++)
	{
	  rA[i] = a[iterf2c (i, sa[0], sa[1])];

	  /*Invalid Bias operation control */

	  if (rinf != 0 || rsup != 0)
	    {
	      if (rA[i] == 0 || fabs (rA[i]) == HUGE_VAL)
		{
		  rinf++;
		  rsup++;
		}
	    }
	}

      if (rinf >= 1 && rsup >= 1)

	{
	  /*Case with invalid Bias operations */

	  DIYMulMRMI (pC, rA, pB, sa, sb);

	  for (i = 0; i < sa[0] * sb[1]; i++)
	    {
	      cinf[iterf2c (i, sa[0], sb[1])] = pC[i].inf;
	      csup[iterf2c (i, sa[0], sb[1])] = pC[i].sup;
	    }

	  free (pB);
	  free (pC);
	  free (rA);
	  *flag = 2;
	  return 0;
	}

      else
	{
	  /*General valid case */

	  BiasMulMRMI (pC, rA, pB, sa[0], sa[1], sb[1]);

	  for (i = 0; i < sa[0] * sb[1]; i++)
	    {
	      cinf[iterf2c (i, sa[0], sb[1])] = pC[i].inf;
	      csup[iterf2c (i, sa[0], sb[1])] = pC[i].sup;
	    }

	  free (pB);
	  free (pC);
	  free (rA);
	  *flag = 2;
	  return 0;
	}
    }
}


int
biasRmulII4sci (const double *ainf, const double *asup, const int *sa,
		const double *binf, const double *bsup, const int *sb,
		double *cinf, double *csup, int *flag)
{
  //feclearexcept (FE_ALL_EXCEPT);

  int i, t[1];
  BIASINTERVAL *pA;
  BIASINTERVAL *pB;
  BIASINTERVAL *pC;

  /*Invalid Bias operation control constants */
  int rinf = 0;
  int rsup = 0;

  double ai[1];
  double as[1];
  double bi[1];
  double bs[1];

  /* Case with different dimensions */
  if (sa[0] != sb[0] || sa[1] != sb[1])
    {
      /*mul a single interval */

      if (sb[0] == 1 && sb[1] == 1)
	{
	  pA =
	    (BIASINTERVAL *) malloc (sa[0] * sa[1] * sizeof (BIASINTERVAL));
	  pB =
	    (BIASINTERVAL *) malloc (sb[0] * sb[1] * sizeof (BIASINTERVAL));
	  pC =
	    (BIASINTERVAL *) malloc (max (sa[0] * sa[1], sb[0] * sb[1]) *
				     sizeof (BIASINTERVAL));
	  Check (pA, pB, pC);
	  Check (pB, pA, pC);
	  Check (pC, pA, pB);

	  pB[0].inf = binf[0];
	  pB[0].sup = bsup[0];

	  /*Invalid Bias operation control */

	  if (pB[0].inf == 0 || fabs (pB[0].inf) == HUGE_VAL)
	    {
	      rinf++;
	    }

	  if (pB[0].sup == 0 || fabs (pB[0].sup) == HUGE_VAL)
	    {
	      rsup++;
	    }

	  for (i = 0; i < sa[0] * sa[1]; i++)
	    {
	      pA[i].inf = ainf[iterf2c (i, sa[0], sa[1])];
	      pA[i].sup = asup[iterf2c (i, sa[0], sa[1])];

	      /*Invalid Bias operation control */

	      if (rinf != 0 || rsup != 0)
		{
		  if (pA[i].inf == 0 || fabs (pA[i].inf) == HUGE_VAL)
		    {
		      rinf++;
		    }

		  if (pA[i].sup == 0 || fabs (pA[i].sup) == HUGE_VAL)
		    {
		      rsup++;
		    }
		}
	    }

	  if (rinf >= 1 && rsup >= 1)

	    {
	      /*Case with invalid Bias operations */

	      DIYMulIMI (pC, &(pB[0]), pA, sb);

	      for (i = 0; i < sa[0] * sa[1]; i++)
		{
		  cinf[iterf2c (i, sa[0], sa[1])] = pC[i].inf;
		  csup[iterf2c (i, sa[0], sa[1])] = pC[i].sup;

		}

	      free (pA);
	      free (pB);
	      free (pC);
	      *flag = 1;
	      return 1;
	    }

	  else
	    {
	      /*General valid case */

	      BiasMulIMI (pC, &(pB[0]), pA, sa[0], sa[1]);

	      for (i = 0; i < sa[0] * sa[1]; i++)
		{
		  cinf[iterf2c (i, sa[0], sa[1])] = pC[i].inf;
		  csup[iterf2c (i, sa[0], sa[1])] = pC[i].sup;
		}

	      free (pA);
	      free (pB);
	      free (pC);
	      *flag = 1;
	      return 1;
	    }
	}

      else if (sa[0] == 1 && sa[1] == 1)
	{
	  pA =
	    (BIASINTERVAL *) malloc (sa[0] * sa[1] * sizeof (BIASINTERVAL));
	  pB =
	    (BIASINTERVAL *) malloc (sb[0] * sb[1] * sizeof (BIASINTERVAL));
	  pC =
	    (BIASINTERVAL *) malloc (max (sa[0] * sa[1], sb[0] * sb[1]) *
				     sizeof (BIASINTERVAL));
	  Check (pA, pB, pC);
	  Check (pB, pA, pC);
	  Check (pC, pA, pB);

	  pA[0].inf = ainf[0];
	  pA[0].sup = asup[0];

	  /*Invalid Bias operation control */

	  if (pA[0].inf == 0 || fabs (pA[0].inf) == HUGE_VAL)
	    {
	      rinf++;
	    }

	  if (pA[0].sup == 0 || fabs (pA[0].sup) == HUGE_VAL)
	    {
	      rsup++;
	    }

	  for (i = 0; i < sb[0] * sb[1]; i++)
	    {
	      pB[i].inf = binf[iterf2c (i, sb[0], sb[1])];
	      pB[i].sup = bsup[iterf2c (i, sb[0], sb[1])];

	      /*Invalid Bias operation control */

	      if (rinf != 0 || rsup != 0)
		{
		  if (pB[i].inf == 0 || fabs (pB[i].inf) == HUGE_VAL)
		    {
		      rinf++;
		    }

		  if (pB[i].sup == 0 || fabs (pB[i].sup) == HUGE_VAL)
		    {
		      rsup++;
		    }
		}
	    }

	  if (rinf >= 1 && rsup >= 1)
	    {
	      /*Case with invalid Bias operations */

	      DIYMulIMI (pC, &(pA[0]), pB, sa);

	      for (i = 0; i < sa[0] * sa[1]; i++)
		{
		  cinf[iterf2c (i, sa[0], sa[1])] = pC[i].inf;
		  csup[iterf2c (i, sa[0], sa[1])] = pC[i].sup;
		}

	      free (pA);
	      free (pB);
	      free (pC);
	      *flag = 1;
	      return 1;
	    }

	  else
	    {
	      /*General valid case */

	      BiasMulIMI (pC, &(pA[0]), pB, sb[0], sb[1]);

	      for (i = 0; i < sb[0] * sb[1]; i++)
		{
		  cinf[iterf2c (i, sb[0], sb[1])] = pC[i].inf;
		  csup[iterf2c (i, sb[0], sb[1])] = pC[i].sup;
		}

	      free (pA);
	      free (pB);
	      free (pC);
	      *flag = 1;
	      return 1;
	    }
	}

      /* Impossible operation */

      else
	{
	  //printf (" INTERVAL Dimension error between matrix");
	  *flag = -1;
	  return -1;
	}
    }

  /*Case with similar dimensions */
  else
    {
      pA = (BIASINTERVAL *) malloc (sa[0] * sa[1] * sizeof (BIASINTERVAL));
      pB = (BIASINTERVAL *) malloc (sb[0] * sb[1] * sizeof (BIASINTERVAL));
      pC =
	(BIASINTERVAL *) malloc (max (sa[0] * sa[1], sb[0] * sb[1]) *
				 sizeof (BIASINTERVAL));
      Check (pA, pB, pC);
      Check (pB, pA, pC);
      Check (pC, pA, pB);

      for (i = 0; i < sa[0] * sa[1]; i++)
	{
	  /*Invalid Bias operation control */

	  *ai = ainf[iterf2c (i, sb[0], sb[1])];
	  *as = asup[iterf2c (i, sb[0], sb[1])];
	  *bi = binf[iterf2c (i, sb[0], sb[1])];
	  *bs = bsup[iterf2c (i, sb[0], sb[1])];

	  checkmul4sci (ai, as, bi, bs, t);

	  if (*t != 0)
	    {
	      checkmul4sci (bi, bs, ai, as, t);
	    }

	  if (*t != 0)
	    {
	      pA[i].inf = ainf[iterf2c (i, sb[0], sb[1])];
	      pA[i].sup = asup[iterf2c (i, sb[0], sb[1])];

	      pB[i].inf = binf[iterf2c (i, sb[0], sb[1])];
	      pB[i].sup = bsup[iterf2c (i, sb[0], sb[1])];
	    }
	  else
	    {
	      pA[i].inf = *ai;
	      pA[i].sup = *as;
	      pB[i].inf = *bi;
	      pB[i].sup = *bs;
	    }

	  BiasMulII (&(pC[i]), &(pA[i]), &(pB[i]));

	  cinf[iterf2c (i, sb[0], sb[1])] = pC[i].inf;
	  csup[iterf2c (i, sb[0], sb[1])] = pC[i].sup;
	}

      free (pA);
      free (pB);
      free (pC);
      *flag = 0;
      return 0;
    }
}



int
biasRmulIR4sci (const double *ainf, const double *asup, const int *sa,
		const double *b, const int *sb,
		double *cinf, double *csup, int *flag)
{
  //feclearexcept (FE_ALL_EXCEPT);

  int i, t[1];
  double ai[1], bi[1], as[1], bs[1];

  BIASINTERVAL pB;
  BIASINTERVAL *pA;
  BIASINTERVAL *pC;
  double *rB;

  /*Invalid Bias operation control constants */
  int rinf = 0;
  int rsup = 0;

  /* Case with different dimensions */
  if (sa[0] != sb[0] || sa[1] != sb[1])
    {
      /*Rmul a single real */

      if (sb[0] == 1 && sb[1] == 1)
	{
	  /*Invalid Bias operation control */

	  if (b[0] == 0 || fabs (b[0]) == HUGE_VAL)
	    {
	      rinf++;
	      rsup++;
	    }

	  pA =
	    (BIASINTERVAL *) malloc (sa[0] * sa[1] * sizeof (BIASINTERVAL));
	  pC =
	    (BIASINTERVAL *) malloc (max (sa[0] * sa[1], sb[0] * sb[1]) *
				     sizeof (BIASINTERVAL));
	  See (pA, pC);
	  See (pC, pA);

	  for (i = 0; i < sa[0] * sa[1]; i++)
	    {
	      pA[i].inf = ainf[iterf2c (i, sa[0], sa[1])];
	      pA[i].sup = asup[iterf2c (i, sa[0], sa[1])];

	      /*Invalid Bias operation control */

	      if (rinf != 0 || rsup != 0)
		{
		  if (pA[i].inf == 0 || fabs (pA[i].inf) == HUGE_VAL)
		    {
		      rinf++;
		    }

		  if (pA[i].sup == 0 || fabs (pA[i].sup) == HUGE_VAL)
		    {
		      rsup++;
		    }
		}
	    }

	  if (rinf > 1 && rsup > 1)
	    {
	      /*Case with invalid Bias operations */

	      DIYMulRMI (pC, b, pA, sa);

	      for (i = 0; i < sa[0] * sa[1]; i++)
		{
		  cinf[iterf2c (i, sa[0], sa[1])] = pC[i].inf;
		  csup[iterf2c (i, sa[0], sa[1])] = pC[i].sup;
		}

	      free (pA);
	      free (pC);
	      *flag = 1;
	      return 1;
	    }

	  else
	    {
	      /*General valid case */

	      BiasMulRMI (pC, &(b[0]), pA, sa[0], sa[1]);

	      for (i = 0; i < sa[0] * sa[1]; i++)
		{
		  cinf[iterf2c (i, sa[0], sa[1])] = pC[i].inf;
		  csup[iterf2c (i, sa[0], sa[1])] = pC[i].sup;
		}

	      free (pA);
	      free (pC);
	      *flag = 1;
	      return 1;
	    }
	}

      /*Rmul a single interval */

      else if (sa[0] == 1 && sa[1] == 1)
	{
	  pA =
	    (BIASINTERVAL *) malloc (sa[0] * sa[1] * sizeof (BIASINTERVAL));
	  rB = (double *) malloc (sb[0] * sb[1] * sizeof (double));
	  pC =
	    (BIASINTERVAL *) malloc (max (sa[0] * sa[1], sb[0] * sb[1]) *
				     sizeof (BIASINTERVAL));
	  Check (pA, rB, pC);
	  Check (rB, pA, pC);
	  Check (pC, pA, rB);

	  pA[0].inf = ainf[0];
	  pA[0].sup = asup[0];

	  /*Invalid Bias operation control */

	  if (pA[0].inf == 0 || fabs (pA[0].inf) == HUGE_VAL)
	    {
	      rinf++;
	    }

	  if (pA[0].sup == 0 || fabs (pA[0].sup) == HUGE_VAL)
	    {
	      rsup++;
	    }

	  for (i = 0; i < sb[0] * sb[1]; i++)
	    {
	      rB[i] = b[iterf2c (i, sb[0], sb[1])];

	      if (rinf != 0 || rsup != 0)
		{
		  if (b[i] == 0 || fabs (b[i]) == HUGE_VAL)
		    {
		      rinf++;
		      rsup++;
		    }
		}
	    }

	  if (rinf >= 2 && rsup >= 2)
	    {
	      /*Case with invalid Bias operations */

	      DIYMulIMR (pC, &(pA[0]), rB, sb);

	      for (i = 0; i < sb[0] * sb[1]; i++)
		{
		  cinf[iterf2c (i, sb[0], sb[1])] = pC[i].inf;
		  csup[iterf2c (i, sb[0], sb[1])] = pC[i].sup;
		}

	      free (pA);
	      free (pC);
	      free (rB);
	      *flag = 1;
	      return 1;
	    }

	  else
	    {
	      /*General valid case */

	      BiasMulIMR (pC, &(pA[0]), rB, sb[0], sb[1]);

	      for (i = 0; i < sb[0] * sb[1]; i++)
		{
		  cinf[iterf2c (i, sb[0], sb[1])] = pC[i].inf;
		  csup[iterf2c (i, sb[0], sb[1])] = pC[i].sup;
		}

	      free (pA);
	      free (pC);
	      *flag = 1;
	      return 1;
	    }
	}

      /*Impossible operation */
      else
	{
	  //printf (" INTERVAL Dimension error between matrix");
	  *flag = -1;
	  return -1;
	}
    }

  /*Case with similar dimensions */
  else
    {
      pA = (BIASINTERVAL *) malloc (sa[0] * sa[1] * sizeof (BIASINTERVAL));
      pC =
	(BIASINTERVAL *) malloc (max (sa[0] * sa[1], sb[0] * sb[1]) *
				 sizeof (BIASINTERVAL));
      See (pA, pC);
      See (pC, pA);

      for (i = 0; i < sa[0] * sa[1]; i++)
	{
	  /*Invalid Bias operation control */

	  *ai = ainf[iterf2c (i, sb[0], sb[1])];
	  *as = asup[iterf2c (i, sb[0], sb[1])];
	  *bi = b[iterf2c (i, sb[0], sb[1])];
	  *bs = b[iterf2c (i, sb[0], sb[1])];

	  checkmul4sci (ai, as, bi, bs, t);

	  if (*t != 0)
	    {
	      checkmul4sci (bi, bs, ai, as, t);
	    }

	  if (*t != 0)
	    {
	      pA[i].inf = ainf[iterf2c (i, sb[0], sb[1])];
	      pA[i].sup = asup[iterf2c (i, sb[0], sb[1])];
	      BiasMulIR (&(pC[i]), &(pA[i]), &(b[iterf2c (i, sb[0], sb[1])]));
	    }
	  else
	    {
	      pA[i].inf = *ai;
	      pA[i].sup = *as;
	      pB.inf = *bi;
	      pB.sup = *bs;
	      BiasMulII (&(pC[i]), &(pA[i]), &(pB));
	    }

	  cinf[iterf2c (i, sa[0], sa[1])] = pC[i].inf;
	  csup[iterf2c (i, sa[0], sa[1])] = pC[i].sup;
	}

      free (pA);
      free (pC);
      *flag = 0;
      return 0;
    }
}


int
biasRLmulII4sci (const double *ainf, const double *asup, const int *sa,
		 const double *binf, const double *bsup, const int *sb,
		 double *cinf, double *csup, int *flag)
{

  //feclearexcept (FE_ALL_EXCEPT);

  int i, j;
  BIASINTERVAL *pA;
  BIASINTERVAL *pB;
  BIASINTERVAL *pC;
  BIASINTERVAL *pE;

  /*Invalid Bias operation control constants */
  int rinf = 0;
  int rsup = 0;

  pA = (BIASINTERVAL *) malloc (sa[0] * sa[1] * sizeof (BIASINTERVAL));
  pB = (BIASINTERVAL *) malloc (sb[0] * sb[1] * sizeof (BIASINTERVAL));
  pC =
    (BIASINTERVAL *) malloc (sa[0] * sa[1] * sb[0] * sb[1] *
			     sizeof (BIASINTERVAL));
  Check (pA, pB, pC);
  Check (pB, pA, pC);
  Check (pC, pA, pB);

  for (i = 0; i < sa[0] * sa[1]; i++)
    {
      pA[i].inf = ainf[iterf2c (i, sa[0], sa[1])];
      pA[i].sup = asup[iterf2c (i, sa[0], sa[1])];

      /*Invalid Bias operation control */

      if (pA[i].inf == 0 || fabs (pA[i].inf) == HUGE_VAL)
	{
	  rinf++;
	}

      if (pA[i].sup == 0 || fabs (pA[i].sup) == HUGE_VAL)
	{
	  rsup++;
	}
    }

  for (i = 0; i < sb[0] * sb[1]; i++)
    {
      pB[i].inf = binf[iterf2c (i, sb[0], sb[1])];
      pB[i].sup = bsup[iterf2c (i, sb[0], sb[1])];

      if (rinf != 0 || rsup != 0)
	{
	  if (pB[i].inf == 0 || fabs (pB[i].inf) == HUGE_VAL)
	    {
	      rinf++;
	    }

	  if (pB[i].sup == 0 || fabs (pB[i].sup) == HUGE_VAL)
	    {
	      rsup++;
	    }
	}
    }

  if (rinf >= 1 && rsup >= 1)
    {
      /*Case with invalid Bias operations */

      DIYKMulMIMI (pC, pA, pB, sa, sb);

      for (i = 0; i < sa[0] * sa[1] * sb[0] * sb[1]; i++)
	{

	  cinf[iterf2c (i, sa[0] * sb[0], sa[1] * sb[1])] = pC[i].inf;
	  csup[iterf2c (i, sa[0] * sb[0], sa[1] * sb[1])] = pC[i].sup;
	}

      free (pA);
      free (pC);
      free (pB);
      *flag = 2;
      return 1;
    }

  else
    {
      /*General valid case */

      for (i = 0; i < sa[0] * sa[1]; i++)
	{
	  pE =
	    (BIASINTERVAL *) malloc (sb[0] * sb[1] * sizeof (BIASINTERVAL));

	  BiasMulIMI (pE, &(pA[i]), pB, sb[0], sb[1]);

	  for (j = 0; j < sb[0] * sb[1]; j++)
	    {
	      pC[j % sb[1] + (i % sa[1]) * sb[1]
		 + (j / sb[1] + (i / sa[1]) * sb[0]) * sa[1] * sb[1]] = pE[j];
	    }
	  free (pE);
	}

      for (i = 0; i < sa[0] * sa[1] * sb[0] * sb[1]; i++)
	{

	  cinf[iterf2c (i, sa[0] * sb[0], sa[1] * sb[1])] = pC[i].inf;
	  csup[iterf2c (i, sa[0] * sb[0], sa[1] * sb[1])] = pC[i].sup;
	}

      free (pA);
      free (pC);
      free (pB);
      *flag = 2;
      return 1;
    }
}


int
biasRLmulIR4sci (const double *ainf, const double *asup, const int *sa,
		 const double *b, const int *sb,
		 double *cinf, double *csup, int *flag)
{

  //feclearexcept (FE_ALL_EXCEPT);

  int i, j;
  BIASINTERVAL *pA;
  double *rB;
  BIASINTERVAL *pC;
  BIASINTERVAL *pE;

  /*Invalid Bias operation control constants */
  int rinf = 0;
  int rsup = 0;

  pA = (BIASINTERVAL *) malloc (sa[0] * sa[1] * sizeof (BIASINTERVAL));
  rB = (double *) malloc (sb[0] * sb[1] * sizeof (double));
  pC =
    (BIASINTERVAL *) malloc (sa[0] * sa[1] * sb[0] * sb[1] *
			     sizeof (BIASINTERVAL));
  Check (pA, rB, pC);
  Check (rB, pA, pC);
  Check (pC, pA, rB);

  for (i = 0; i < sa[0] * sa[1]; i++)
    {
      pA[i].inf = ainf[iterf2c (i, sa[0], sa[1])];
      pA[i].sup = asup[iterf2c (i, sa[0], sa[1])];

      /*Invalid Bias operation control */

      if (pA[i].inf == 0 || fabs (pA[i].inf) == HUGE_VAL)
	{
	  rinf++;
	}

      if (pA[i].sup == 0 || fabs (pA[i].sup) == HUGE_VAL)
	{
	  rsup++;
	}
    }

  for (i = 0; i < sb[0] * sb[1]; i++)
    {
      rB[i] = b[iterf2c (i, sb[0], sb[1])];

      if (rinf != 0 || rsup != 0)
	{
	  if (b[i] == 0 || fabs (b[i]) == HUGE_VAL)
	    {
	      rinf++;
	      rsup++;
	    }
	}
    }

  if (rinf >= 2 && rsup >= 2)
    {
      /*Case with invalid Bias operations */

      DIYKMulMIMR (pC, pA, rB, sa, sb);

      for (i = 0; i < sa[0] * sa[1] * sb[0] * sb[1]; i++)
	{
	  cinf[iterf2c (i, sa[0] * sb[0], sa[1] * sb[1])] = pC[i].inf;
	  csup[iterf2c (i, sa[0] * sb[0], sa[1] * sb[1])] = pC[i].sup;
	}

      free (pA);
      free (pC);
      free (rB);
      *flag = 2;
      return 1;
    }

  else
    {
      /*General valid case */

      for (i = 0; i < sa[0] * sa[1]; i++)
	{
	  pE =
	    (BIASINTERVAL *) malloc (sb[0] * sb[1] * sizeof (BIASINTERVAL));

	  BiasMulIMR (pE, &(pA[i]), rB, sb[0], sb[1]);

	  for (j = 0; j < sb[0] * sb[1]; j++)
	    {
	      pC[j % sb[1] + (i % sa[1]) * sb[1]
		 + (j / sb[1] + (i / sa[1]) * sb[0]) * sa[1] * sb[1]] = pE[j];
	    }
	  free (pE);
	}

      for (i = 0; i < sa[0] * sa[1] * sb[0] * sb[1]; i++)
	{
	  cinf[iterf2c (i, sa[0] * sb[0], sa[1] * sb[1])] = pC[i].inf;
	  csup[iterf2c (i, sa[0] * sb[0], sa[1] * sb[1])] = pC[i].sup;
	}

      free (pA);
      free (pC);
      free (rB);
      *flag = 2;
      return 1;
    }
}


int
biasRLmulRI4sci (const double *a, const int *sa,
		 const double *binf, const double *bsup, const int *sb,
		 double *cinf, double *csup, int *flag)
{
  //feclearexcept (FE_ALL_EXCEPT);

  int i, j;
  double *rA;
  BIASINTERVAL *pB;
  BIASINTERVAL *pC;
  BIASINTERVAL *pE;

  /*Invalid Bias operation control constants */
  int rinf = 0;
  int rsup = 0;

  rA = (double *) malloc (sa[0] * sa[1] * sizeof (double));
  pB = (BIASINTERVAL *) malloc (sb[0] * sb[1] * sizeof (BIASINTERVAL));
  pC =
    (BIASINTERVAL *) malloc (sa[0] * sa[1] * sb[0] * sb[1] *
			     sizeof (BIASINTERVAL));
  Check (rA, pB, pC);
  Check (pB, rA, pC);
  Check (pC, rA, pB);

  for (i = 0; i < sa[0] * sa[1]; i++)
    {
      rA[i] = a[iterf2c (i, sa[0], sa[1])];

      if (a[i] == 0 || fabs (a[i]) == HUGE_VAL)
	{
	  rinf++;
	  rsup++;
	}
    }

  for (i = 0; i < sb[0] * sb[1]; i++)
    {
      pB[i].inf = binf[iterf2c (i, sb[0], sb[1])];
      pB[i].sup = bsup[iterf2c (i, sb[0], sb[1])];

      if (rinf != 0 || rsup != 0)
	{
	  if (pB[i].inf == 0 || fabs (pB[i].inf) == HUGE_VAL)
	    {
	      rinf++;
	    }

	  if (pB[i].sup == 0 || fabs (pB[i].sup) == HUGE_VAL)
	    {
	      rsup++;
	    }
	}
    }

  if (rinf >= 1 && rsup >= 1)
    {
      /*Case with invalid Bias operations */

      DIYKMulMRMI (pC, rA, pB, sa, sb);

      for (i = 0; i < sa[0] * sa[1] * sb[0] * sb[1]; i++)
	{
	  cinf[iterf2c (i, sa[0] * sb[0], sa[1] * sb[1])] = pC[i].inf;
	  csup[iterf2c (i, sa[0] * sb[0], sa[1] * sb[1])] = pC[i].sup;
	}

      free (rA);
      free (pC);
      free (pB);
      *flag = 2;
      return 1;
    }

  else
    {
      /* General valid case */

      for (i = 0; i < sa[0] * sa[1]; i++)
	{
	  pE =
	    (BIASINTERVAL *) malloc (sb[0] * sb[1] * sizeof (BIASINTERVAL));

	  BiasMulRMI (pE, &(rA[i]), pB, sb[0], sb[1]);

	  for (j = 0; j < sb[0] * sb[1]; j++)
	    {
	      pC[j % sb[1] + (i % sa[1]) * sb[1]
		 + (j / sb[1] + (i / sa[1]) * sb[0]) * sa[1] * sb[1]] = pE[j];

	    }
	  free (pE);
	}

      for (i = 0; i < sa[0] * sa[1] * sb[0] * sb[1]; i++)
	{
	  cinf[iterf2c (i, sa[0] * sb[0], sa[1] * sb[1])] = pC[i].inf;
	  csup[iterf2c (i, sa[0] * sb[0], sa[1] * sb[1])] = pC[i].sup;
	}

      free (rA);
      free (pC);
      free (pB);
      *flag = 2;
      return 1;
    }
}

int
biasdivII4sci (const double *ainf, const double *asup, const int *sa,
	       const double *binf, const double *bsup, const int *sb,
	       double *cinf, double *csup, int *flag)
{
  //feclearexcept (FE_ALL_EXCEPT);

  int i;
  BIASINTERVAL *pA;
  BIASINTERVAL *pB;
  BIASINTERVAL *pC;

  /*Invalid Bias Operation control constants */
  int r0 = 0;
  int rinf = 0;
  int rsup = 0;


  /* Case with different dimensions */
  if (sa[0] != sb[0] || sa[1] != sb[1])
    {
      /*div by a single interval */

      if (sb[0] == 1 && sb[1] == 1)
	{
	  pA =
	    (BIASINTERVAL *) malloc (sa[0] * sa[1] * sizeof (BIASINTERVAL));
	  pB =
	    (BIASINTERVAL *) malloc (sb[0] * sb[1] * sizeof (BIASINTERVAL));
	  pC =
	    (BIASINTERVAL *) malloc (max (sa[0] * sa[1], sb[0] * sb[1]) *
				     sizeof (BIASINTERVAL));
	  Check (pA, pB, pC);
	  Check (pB, pA, pC);
	  Check (pC, pA, pB);

	  pB[0].inf = binf[0];
	  pB[0].sup = bsup[0];

	  /*Invalid Bias operation control */
	  if (pB[0].inf < 0 && pB[0].sup > 0)
	    {
	      r0++;
	    }

	  if (pB[0].inf == 0 || fabs (pB[0].inf) == HUGE_VAL)
	    {
	      rinf++;
	    }

	  if (pB[0].sup == 0 || fabs (pB[0].sup) == HUGE_VAL)
	    {
	      rsup++;
	    }

	  for (i = 0; i < sa[0] * sa[1]; i++)
	    {
	      pA[i].inf = ainf[iterf2c (i, sa[0], sa[1])];
	      pA[i].sup = asup[iterf2c (i, sa[0], sa[1])];
	    }

	  /* Division by zero */

	  if (r0 != 0)
	    {
	      for (i = 0; i < sa[0] * sa[1]; i++)
		{
		  cinf[iterf2c (i, sa[0], sa[1])] = -HUGE_VAL;
		  csup[iterf2c (i, sa[0], sa[1])] = HUGE_VAL;
		}

	      free (pA);
	      free (pB);
	      free (pC);
	      *flag = 1;
	      return 1;
	    }

	  /* Invalid Bias Operation Control */

	  if (rinf != 0 || rsup != 0)
	    {
	      for (i = 0; i < sa[0] * sa[1]; i++)
		{
		  DIYDivII (&(pC[i]), &(pA[i]), &(pB[0]));
		}
	    }
	  else
	    {
	      /*General valid case */
	      BiasDivMII (pC, pA, &(pB[0]), sa[0], sa[1]);
	    }

	  for (i = 0; i < sa[0] * sa[1]; i++)
	    {
	      cinf[iterf2c (i, sa[0], sa[1])] = pC[i].inf;
	      csup[iterf2c (i, sa[0], sa[1])] = pC[i].sup;
	    }

	  free (pA);
	  free (pB);
	  free (pC);
	  *flag = 1;
	  return 1;

	}
      else if (sa[0] == 1 && sa[1] == 1)
	{
	  pA =
	    (BIASINTERVAL *) malloc (sa[0] * sa[1] * sizeof (BIASINTERVAL));
	  pB =
	    (BIASINTERVAL *) malloc (sb[0] * sb[1] * sizeof (BIASINTERVAL));
	  pC =
	    (BIASINTERVAL *) malloc (max (sa[0] * sa[1], sb[0] * sb[1]) *
				     sizeof (BIASINTERVAL));
	  Check (pA, pB, pC);
	  Check (pB, pA, pC);
	  Check (pC, pA, pB);

	  pA[0].inf = ainf[0];
	  pA[0].sup = asup[0];

	  for (i = 0; i < sb[0] * sb[1]; i++)
	    {
	      pB[i].inf = binf[iterf2c (i, sb[0], sb[1])];
	      pB[i].sup = bsup[iterf2c (i, sb[0], sb[1])];

	      /*Divsion by zero */
	      if (pB[i].inf < 0 && pB[i].sup > 0)
		{
		  pC[i].inf = -HUGE_VAL;
		  pC[i].sup = HUGE_VAL;
		}

	      /* Invalid Bias Operation Control */

	      else if ((pB[i].inf == 0 || fabs (pB[i].inf) == HUGE_VAL) ||
		       (pB[i].sup == 0 || fabs (pB[i].sup) == HUGE_VAL))
		{
		  DIYDivII (&(pC[i]), &(pA[0]), &(pB[i]));
		}

	      else
		{
		  /*General Valid Case */
		  BiasDivII (&(pC[i]), &(pA[0]), &(pB[i]));
		}

	      cinf[iterf2c (i, sb[0], sb[1])] = pC[i].inf;
	      csup[iterf2c (i, sb[0], sb[1])] = pC[i].sup;
	    }

	  free (pA);
	  free (pB);
	  free (pC);
	  *flag = 1;
	  return 1;
	}

      /* Impossible operation */

      else
	{
	  //printf (" INTERVAL Dimension error between matrix");
	  *flag = -1;
	  return -1;
	}
    }

  /*Case with similar dimensions */
  else
    {
      if (sa[0] * sa[1] == 1)
	{
	  pA =
	    (BIASINTERVAL *) malloc (sa[0] * sa[1] * sizeof (BIASINTERVAL));
	  pB =
	    (BIASINTERVAL *) malloc (sb[0] * sb[1] * sizeof (BIASINTERVAL));
	  pC =
	    (BIASINTERVAL *) malloc (max (sa[0] * sa[1], sb[0] * sb[1]) *
				     sizeof (BIASINTERVAL));
	  Check (pA, pB, pC);
	  Check (pB, pA, pC);
	  Check (pC, pA, pB);

	  pA[0].inf = ainf[0];
	  pA[0].sup = asup[0];

	  pB[0].inf = binf[0];
	  pB[0].sup = bsup[0];

	  /*Divsion by zero */

	  if (pB[0].inf < 0 && pB[0].sup > 0)
	    {
	      pC[0].inf = -HUGE_VAL;
	      pC[0].sup = HUGE_VAL;
	    }

	  /* Invalid Bias Operation Control */

	  else if ((pB[0].inf == 0 || fabs (pB[0].inf) == HUGE_VAL) ||
		   (pB[0].sup == 0 || fabs (pB[0].sup) == HUGE_VAL))
	    {
	      DIYDivII (&(pC[0]), &(pA[0]), &(pB[0]));
	    }

	  else
	    {
	      /*General valid Case */
	      BiasDivII (&(pC[0]), &(pA[0]), &(pB[0]));
	    }

	  cinf[0] = pC[0].inf;
	  csup[0] = pC[0].sup;

	  free (pA);
	  free (pB);
	  free (pC);
	  *flag = 0;
	  return 0;
	}

      else
	{
	  // Not implemented yet !!
	  *flag = -3;
	  return 0;
	}
    }
}



int
biasdivIR4sci (const double *ainf, const double *asup, const int *sa,
	       const double *b, const int *sb,
	       double *cinf, double *csup, int *flag)
{
  //feclearexcept (FE_ALL_EXCEPT);

  int i;
  BIASINTERVAL *pA;
  BIASINTERVAL *pC;
  double *rB;

  /* Case with different dimensions */
  if (sa[0] != sb[0] || sa[1] != sb[1])
    {
      /*div by a single real */

      if (sb[0] == 1 && sb[1] == 1)
	{
	  pA =
	    (BIASINTERVAL *) malloc (sa[0] * sa[1] * sizeof (BIASINTERVAL));
	  pC =
	    (BIASINTERVAL *) malloc (max (sa[0] * sa[1], sb[0] * sb[1]) *
				     sizeof (BIASINTERVAL));
	  See (pA, pC);
	  See (pC, pA);

	  for (i = 0; i < sa[0] * sa[1]; i++)
	    {
	      pA[i].inf = ainf[iterf2c (i, sa[0], sa[1])];
	      pA[i].sup = asup[iterf2c (i, sa[0], sa[1])];
	    }

	  /*Division by zero */

	  if (b[0] == 0)
	    {
	      for (i = 0; i < sa[0] * sa[1]; i++)
		{
		  pC[i].inf = -HUGE_VAL;
		  pC[i].sup = HUGE_VAL;
		}
	    }

	  /*Invalid Bias Operations Control */

	  else if (fabs (b[0]) == HUGE_VAL)
	    {
	      for (i = 0; i < sa[0] * sa[1]; i++)
		{
		  pC[i].inf = 0;
		  pC[i].sup = 0;
		}
	    }

	  else
	    /*General Valid case */
	    {
	      BiasDivMIR (pC, pA, &(b[0]), sa[0], sa[1]);
	    }

	  for (i = 0; i < sa[0] * sa[1]; i++)
	    {
	      cinf[iterf2c (i, sa[0], sa[1])] = pC[i].inf;
	      csup[iterf2c (i, sa[0], sa[1])] = pC[i].sup;
	    }

	  free (pA);
	  free (pC);
	  *flag = 1;
	  return 1;
	}

      /*Div a single interval */

      else if (sa[0] == 1 && sa[1] == 1)
	{
	  pA =
	    (BIASINTERVAL *) malloc (sa[0] * sa[1] * sizeof (BIASINTERVAL));
	  rB = (double *) malloc (sb[0] * sb[1] * sizeof (double));
	  pC =
	    (BIASINTERVAL *) malloc (max (sa[0] * sa[1], sb[0] * sb[1]) *
				     sizeof (BIASINTERVAL));
	  Check (pA, rB, pC);
	  Check (rB, pA, pC);
	  Check (pC, pA, rB);

	  pA[0].inf = ainf[0];
	  pA[0].sup = asup[0];

	  for (i = 0; i < sb[0] * sb[1]; i++)
	    {
	      rB[i] = b[iterf2c (i, sb[0], sb[1])];

	      /*Division by zero */

	      if (rB[i] == 0)
		{
		  pC[i].inf = -HUGE_VAL;
		  pC[i].sup = HUGE_VAL;
		}

	      /*Invalid Bias operations Control */

	      else if (fabs (rB[i]) == HUGE_VAL)
		{
		  pC[i].inf = 0;
		  pC[i].sup = 0;
		}

	      else
		{
		  /*General Valid Case */
		  BiasDivIR (&(pC[i]), &(pA[0]), &(rB[i]));
		}

	      cinf[iterf2c (i, sb[0], sb[1])] = pC[i].inf;
	      csup[iterf2c (i, sb[0], sb[1])] = pC[i].sup;
	    }

	  free (pA);
	  free (rB);
	  free (pC);
	  *flag = 1;
	  return 1;
	}

      /*Impossible operation */
      else
	{
	  //printf (" INTERVAL Dimension error between matrix");
	  *flag = -1;
	  return -1;
	}
    }

  /*Case with similar dimensions */
  else
    {
      if (sa[0] * sa[1] == 1)
	{
	  pA =
	    (BIASINTERVAL *) malloc (sa[0] * sa[1] * sizeof (BIASINTERVAL));
	  rB = (double *) malloc (sb[0] * sb[1] * sizeof (double));
	  pC =
	    (BIASINTERVAL *) malloc (max (sa[0] * sa[1], sb[0] * sb[1]) *
				     sizeof (BIASINTERVAL));
	  Check (pA, rB, pC);
	  Check (rB, pA, pC);
	  Check (pC, pA, rB);

	  pA[0].inf = ainf[0];
	  pA[0].sup = asup[0];

	  rB[0] = b[0];

	  /*Division by zero */
	  if (b[0] == 0)
	    {
	      pC[0].inf = -HUGE_VAL;
	      pC[0].sup = HUGE_VAL;
	    }

	  /*Invalid Bias Operation Control */
	  else if (fabs (b[0]) == HUGE_VAL)
	    {
	      pC[0].inf = 0;
	      pC[0].sup = 0;
	    }

	  /*General valid Case */
	  else
	    {
	      BiasDivIR (&(pC[0]), &(pA[0]), &(rB[0]));
	    }

	  cinf[0] = pC[0].inf;
	  csup[0] = pC[0].sup;

	  free (pA);
	  free (rB);
	  free (pC);
	  *flag = 0;
	  return 0;
	}

      else
	{
	  // Not implemented yet !!
	  *flag = -3;
	  return 0;

	}
    }
}


int
biasdivRI4sci (const double *a, const int *sa,
	       const double *binf, const double *bsup, const int *sb,
	       double *cinf, double *csup, int *flag)
{
  //feclearexcept (FE_ALL_EXCEPT);


  int i;
  BIASINTERVAL *pB;
  BIASINTERVAL *pC;
  double *rA;

  /*Invalid Bias Operation control constants */
  int r0 = 0;
  int rinf = 0;
  int rsup = 0;

  /* Case with different dimensions */
  if (sa[0] != sb[0] || sa[1] != sb[1])
    {
      /*div by a single interval */

      if (sb[0] == 1 && sb[1] == 1)
	{
	  rA = (double *) malloc (sa[0] * sa[1] * sizeof (double));
	  pB =
	    (BIASINTERVAL *) malloc (sb[0] * sb[1] * sizeof (BIASINTERVAL));
	  pC =
	    (BIASINTERVAL *) malloc (max (sa[0] * sa[1], sb[0] * sb[1]) *
				     sizeof (BIASINTERVAL));
	  Check (rA, pB, pC);
	  Check (pB, rA, pC);
	  Check (pC, rA, pB);

	  pB[0].inf = binf[0];
	  pB[0].sup = bsup[0];

	  /*Invalid Bias operation control */
	  if (pB[0].inf < 0 && pB[0].sup > 0)
	    {
	      r0++;
	    }

	  if (pB[0].inf == 0 || fabs (pB[0].inf) == HUGE_VAL)
	    {
	      rinf++;
	    }

	  if (pB[0].sup == 0 || fabs (pB[0].sup) == HUGE_VAL)
	    {
	      rsup++;
	    }

	  for (i = 0; i < sa[0] * sa[1]; i++)
	    {
	      rA[i] = a[iterf2c (i, sa[0], sa[1])];
	    }

	  /* Division by zero */

	  if (r0 != 0)
	    {
	      for (i = 0; i < sa[0] * sa[1]; i++)
		{
		  cinf[iterf2c (i, sa[0], sa[1])] = -HUGE_VAL;
		  csup[iterf2c (i, sa[0], sa[1])] = HUGE_VAL;
		}

	      free (rA);
	      free (pB);
	      free (pC);
	      *flag = 1;
	      return 1;
	    }

	  /* Invalid Bias Operation Control */

	  if (rinf != 0 || rsup != 0)
	    {
	      for (i = 0; i < sa[0] * sa[1]; i++)
		{
		  DIYDivRI (&(pC[i]), &(rA[i]), &(pB[i]));
		}
	    }
	  else
	    {
	      /*General valid case */
	      BiasDivMRI (pC, rA, &(pB[0]), sa[0], sa[1]);
	    }

	  for (i = 0; i < sa[0] * sa[1]; i++)
	    {
	      cinf[iterf2c (i, sa[0], sa[1])] = pC[i].inf;
	      csup[iterf2c (i, sa[0], sa[1])] = pC[i].sup;
	    }

	  free (rA);
	  free (pB);
	  free (pC);
	  *flag = 1;
	  return 1;

	}

      /*Div a single real */

      else if (sa[0] == 1 && sa[1] == 1)
	{
	  rA = (double *) malloc (sa[0] * sa[1] * sizeof (double));
	  pB =
	    (BIASINTERVAL *) malloc (sb[0] * sb[1] * sizeof (BIASINTERVAL));
	  pC =
	    (BIASINTERVAL *) malloc (max (sa[0] * sa[1], sb[0] * sb[1]) *
				     sizeof (BIASINTERVAL));
	  Check (rA, pB, pC);
	  Check (pB, rA, pC);
	  Check (pC, rA, pB);

	  rA[0] = a[0];

	  for (i = 0; i < sb[0] * sb[1]; i++)
	    {
	      pB[i].inf = binf[iterf2c (i, sb[0], sb[1])];
	      pB[i].sup = bsup[iterf2c (i, sb[0], sb[1])];

	      /*Divsion by zero */
	      if (pB[i].inf < 0 && pB[i].sup > 0)
		{
		  pC[i].inf = -HUGE_VAL;
		  pC[i].sup = HUGE_VAL;
		}

	      /* Invalid Bias Operation Control */

	      else if ((pB[i].inf == 0 || fabs (pB[i].inf) == HUGE_VAL) ||
		       (pB[i].sup == 0 || fabs (pB[i].sup) == HUGE_VAL))
		{
		  DIYDivRI (&(pC[i]), &(rA[0]), &(pB[i]));
		}

	      else
		{
		  /*General Valid Case */
		  BiasDivRI (&(pC[i]), &(rA[0]), &(pB[i]));
		}

	      cinf[iterf2c (i, sb[0], sb[1])] = pC[i].inf;
	      csup[iterf2c (i, sb[0], sb[1])] = pC[i].sup;

	    }

	  free (rA);
	  free (pB);
	  free (pC);
	  *flag = 1;
	  return 1;
	}

      /*Impossible operation */
      else
	{
	  //printf (" INTERVAL Dimension error between matrix");
	  *flag = -1;
	  return -1;
	}
    }

  /*Case with similar dimensions */
  else
    {
      if (sa[0] * sa[1] == 1)
	{
	  rA = (double *) malloc (sa[0] * sa[1] * sizeof (double));
	  pB =
	    (BIASINTERVAL *) malloc (sb[0] * sb[1] * sizeof (BIASINTERVAL));
	  pC =
	    (BIASINTERVAL *) malloc (max (sa[0] * sa[1], sb[0] * sb[1]) *
				     sizeof (BIASINTERVAL));
	  Check (rA, pB, pC);
	  Check (pB, rA, pC);
	  Check (pC, rA, pB);

	  rA[0] = a[0];

	  pB[0].inf = binf[0];
	  pB[0].sup = bsup[0];

	  /*Divsion by zero */

	  if (pB[0].inf < 0 && pB[0].sup > 0)
	    {
	      pC[0].inf = -HUGE_VAL;
	      pC[0].sup = HUGE_VAL;
	    }

	  /* Invalid Bias Operation Control */

	  else if ((pB[0].inf == 0 || fabs (pB[0].inf) == HUGE_VAL) ||
		   (pB[0].sup == 0 || fabs (pB[0].sup) == HUGE_VAL))
	    {
	      DIYDivRI (&(pC[0]), &(rA[0]), &(pB[0]));
	    }

	  else
	    {
	      /*General valid Case */
	      BiasDivRI (&(pC[0]), &(rA[0]), &(pB[0]));
	    }

	  cinf[0] = pC[0].inf;
	  csup[0] = pC[0].sup;

	  free (rA);
	  free (pB);
	  free (pC);
	  *flag = 0;
	  return 0;
	}

      else
	{
	  // Not implemented yet !!
	  *flag = -3;
	  return 0;
	}
    }
}


int
biasRdivII4sci (const double *ainf, const double *asup, const int *sa,
		const double *binf, const double *bsup, const int *sb,
		double *cinf, double *csup, int *flag)
{
  //feclearexcept (FE_ALL_EXCEPT);

  int i;
  BIASINTERVAL *pA;
  BIASINTERVAL *pB;
  BIASINTERVAL *pC;

  /* Case with different dimensions */
  if (sa[0] != sb[0] || sa[1] != sb[1])
    {
      biasdivII4sci (ainf, asup, sa, binf, bsup, sb, cinf, csup, flag);
      return 0;
    }

  /*Case with similar dimensions */
  else
    {
      pA = (BIASINTERVAL *) malloc (sa[0] * sa[1] * sizeof (BIASINTERVAL));
      pB = (BIASINTERVAL *) malloc (sb[0] * sb[1] * sizeof (BIASINTERVAL));
      pC =
	(BIASINTERVAL *) malloc (max (sa[0] * sa[1], sb[0] * sb[1]) *
				 sizeof (BIASINTERVAL));
      Check (pA, pB, pC);
      Check (pB, pA, pC);
      Check (pC, pA, pB);

      for (i = 0; i < sa[0] * sa[1]; i++)
	{
	  pA[i].inf = ainf[iterf2c (i, sa[0], sa[1])];
	  pA[i].sup = asup[iterf2c (i, sa[0], sa[1])];

	  pB[i].inf = binf[iterf2c (i, sb[0], sb[1])];
	  pB[i].sup = bsup[iterf2c (i, sb[0], sb[1])];

	  /*Divsion by zero */

	  if (pB[i].inf < 0 && pB[i].sup > 0)
	    {
	      pC[i].inf = -HUGE_VAL;
	      pC[i].sup = HUGE_VAL;
	    }

	  /* Invalid Bias Operation Control */

	  else if ((pB[i].inf == 0 || fabs (pB[i].inf) == HUGE_VAL) ||
		   (pB[i].sup == 0 || fabs (pB[i].sup) == HUGE_VAL))
	    {
	      DIYDivII (&(pC[i]), &(pA[i]), &(pB[i]));
	    }

	  else
	    {
	      /*General valid Case */
	      BiasDivII (&(pC[i]), &(pA[i]), &(pB[i]));
	    }

	  cinf[iterf2c (i, sb[0], sb[1])] = pC[i].inf;
	  csup[iterf2c (i, sb[0], sb[1])] = pC[i].sup;
	}

      free (pA);
      free (pB);
      free (pC);
      *flag = 0;
      return 0;
    }
}



int
biasRdivIR4sci (const double *ainf, const double *asup, const int *sa,
		const double *b, const int *sb,
		double *cinf, double *csup, int *flag)
{
  //feclearexcept (FE_ALL_EXCEPT);

  int i;
  BIASINTERVAL *pA;
  BIASINTERVAL *pC;
  double *rB;

  /* Case with different dimensions */
  if (sa[0] != sb[0] || sa[1] != sb[1])
    {
      biasdivIR4sci (ainf, asup, sa, b, sb, cinf, csup, flag);
      return 0;
    }

  /*Case with similar dimensions */
  else
    {

      pA = (BIASINTERVAL *) malloc (sa[0] * sa[1] * sizeof (BIASINTERVAL));
      rB = (double *) malloc (sb[0] * sb[1] * sizeof (double));
      pC =
	(BIASINTERVAL *) malloc (max (sa[0] * sa[1], sb[0] * sb[1]) *
				 sizeof (BIASINTERVAL));
      Check (pA, rB, pC);
      Check (rB, pA, pC);
      Check (pC, pA, rB);

      for (i = 0; i < sa[0] * sa[1]; i++)
	{
	  pA[i].inf = ainf[iterf2c (i, sa[0], sa[1])];
	  pA[i].sup = asup[iterf2c (i, sa[0], sa[1])];

	  rB[i] = b[iterf2c (i, sb[0], sb[1])];

	  /*Division by zero */
	  if (rB[i] == 0)
	    {
	      pC[i].inf = -HUGE_VAL;
	      pC[i].sup = HUGE_VAL;
	    }

	  /*Invalid Bias Operation Control */
	  else if (fabs (rB[i]) == HUGE_VAL)
	    {
	      pC[i].inf = 0;
	      pC[i].sup = 0;
	    }

	  /*General valid Case */
	  else
	    {
	      BiasDivIR (&(pC[i]), &(pA[i]), &(rB[i]));
	    }

	  cinf[iterf2c (i, sa[0], sa[1])] = pC[i].inf;
	  csup[iterf2c (i, sa[0], sa[1])] = pC[i].sup;
	}

      free (pA);
      free (rB);
      free (pC);
      *flag = 0;
      return 0;
    }
}


int
biasRdivRI4sci (const double *a, const int *sa,
		const double *binf, const double *bsup, const int *sb,
		double *cinf, double *csup, int *flag)
{
  //feclearexcept (FE_ALL_EXCEPT);

  int i;
  BIASINTERVAL *pB;
  BIASINTERVAL *pC;
  double *rA;

  /* Case with different dimensions */
  if (sa[0] != sb[0] || sa[1] != sb[1])
    {
      biasdivRI4sci (a, sa, binf, bsup, sb, cinf, csup, flag);
      return 0;
    }

  /*Case with similar dimensions */
  else
    {
      rA = (double *) malloc (sa[0] * sa[1] * sizeof (double));
      pB = (BIASINTERVAL *) malloc (sb[0] * sb[1] * sizeof (BIASINTERVAL));
      pC =
	(BIASINTERVAL *) malloc (max (sa[0] * sa[1], sb[0] * sb[1]) *
				 sizeof (BIASINTERVAL));
      Check (rA, pB, pC);
      Check (pB, rA, pC);
      Check (pC, rA, pB);

      for (i = 0; i < sa[0] * sa[1]; i++)
	{
	  rA[i] = a[iterf2c (i, sa[0], sa[1])];

	  pB[i].inf = binf[iterf2c (i, sb[0], sb[1])];
	  pB[i].sup = bsup[iterf2c (i, sb[0], sb[1])];

	  /*Divsion by zero */

	  if (pB[i].inf < 0 && pB[i].sup > 0)
	    {
	      pC[i].inf = -HUGE_VAL;
	      pC[i].sup = HUGE_VAL;
	    }

	  /* Invalid Bias Operation Control */

	  else if ((pB[i].inf == 0 || fabs (pB[i].inf) == HUGE_VAL) ||
		   (pB[i].sup == 0 || fabs (pB[i].sup) == HUGE_VAL))
	    {
	      DIYDivRI (&(pC[i]), &(rA[i]), &(pB[i]));
	    }

	  else
	    {
	      /*General valid Case */
	      BiasDivRI (&(pC[i]), &(rA[i]), &(pB[i]));
	    }

	  cinf[iterf2c (i, sb[0], sb[1])] = pC[i].inf;
	  csup[iterf2c (i, sb[0], sb[1])] = pC[i].sup;
	}

      free (pB);
      free (rA);
      free (pC);
      *flag = 0;
      return 0;
    }

}

int
biasrootNI4sci (const double *ainf, const double *asup, const int *sa,
		const int *n, double *cinf, double *csup, int *flag)
{
  //feclearexcept (FE_ALL_EXCEPT);

  int i;
  double zero = 0.0;
  BIASINTERVAL *pA;
  BIASINTERVAL *pC;

  pA = (BIASINTERVAL *) malloc (sa[0] * sa[1] * sizeof (BIASINTERVAL));
  pC = (BIASINTERVAL *) malloc (sa[0] * sa[1] * sizeof (BIASINTERVAL));
  See (pA, pC);
  See (pC, pA);

  BiasFuncInit ();

  for (i = 0; i < sa[0] * sa[1]; i++)
    {
      if (asup[i] < 0)
	{
	  *flag = -1;
	  free (pA);
	  free (pC);
	  return 0;
	}

      pA[i].inf = ainf[i];
      pA[i].sup = asup[i];

      if (BiasInR (&zero, &(pA[i])))
	{
	  pA[i].inf = 0;
	}
    }

  if (n[0] == 2)
    {
      for (i = 0; i < sa[0] * sa[1]; i++)
	{
	  BiasSqrt (&pC[i], &(pA[i]));
	  cinf[i] = pC[i].inf;
	  csup[i] = pC[i].sup;
	}

      *flag = 0;
      free (pA);
      free (pC);
      return 0;
    }

  else
    {
      for (i = 0; i < sa[0] * sa[1]; i++)
	{
	  BiasRoot (&pC[i], &(pA[i]), n[0]);
	  cinf[i] = pC[i].inf;
	  csup[i] = pC[i].sup;
	}
      free (pA);
      free (pC);
      *flag = 0;
      return 0;
    }
}

int
biaspowerNI4sci (const double *ainf, const double *asup, const int *sa,
		 const int *n, double *cinf, double *csup, int *flag)
{
  //feclearexcept (FE_ALL_EXCEPT);

  int i;
  int j;
  int K;
  int nn;
  int Bin[100];
  int Init_check;
  double *dinf;
  double *dsup;
  double *einf;
  double *esup;
  double *finf;
  double *fsup;
  double *ginf;
  double *gsup;

  dinf = (double *) malloc (sa[0] * sa[1] * sizeof (double));
  dsup = (double *) malloc (sa[0] * sa[1] * sizeof (double));
  einf = (double *) malloc (sa[0] * sa[1] * sizeof (double));
  esup = (double *) malloc (sa[0] * sa[1] * sizeof (double));
  finf = (double *) malloc (sa[0] * sa[1] * sizeof (double));
  fsup = (double *) malloc (sa[0] * sa[1] * sizeof (double));
  ginf = (double *) malloc (sa[0] * sa[1] * sizeof (double));
  gsup = (double *) malloc (sa[0] * sa[1] * sizeof (double));
  See (dinf, dsup);
  See (einf, esup);
  See (finf, fsup);
  See (ginf, gsup);

  BiasFuncInit ();
  if (*n == 0)
    {
      *flag = 0;
      free (dinf);
      free (dsup);
      free (einf);
      free (esup);
      free (finf);
      free (fsup);
      free (ginf);
      free (gsup);
      return 0;
    }

  else if (*n == 1)
    {
      for (i = 0; i < sa[0] * sa[1]; i++)
	{
	  cinf[i] = ainf[i];
	  csup[i] = asup[i];
	}
      *flag = 0;
      free (dinf);
      free (dsup);
      free (einf);
      free (esup);
      free (finf);
      free (fsup);
      free (ginf);
      free (gsup);
      return 0;
    }

  else if (*n >= 2)
    {
      // Square matrix
      biasmulII4sci (ainf, asup, sa, ainf, asup, sa, dinf, dsup, flag);

      if (*n == 2)
	{
	  for (i = 0; i < sa[0] * sa[1]; i++)
	    {
	      cinf[i] = dinf[i];
	      csup[i] = dsup[i];
	    }
	  free (dinf);
	  free (dsup);
	  free (einf);
	  free (esup);
	  free (finf);
	  free (fsup);
	  free (ginf);
	  free (gsup);
	  return 0;
	}

      else
	{
	  K = 0;
	  nn = *n;

	  // n written in binary :
	  do
	    {
	      Bin[K] = nn % 2;
	      K++;
	      nn = nn / 2;
	    }
	  while (nn != 0);


	  // n = 3 :
	  if (K == 2)
	    {
	      biasmulII4sci (dinf, dsup, sa, ainf, asup, sa, einf, esup,
			     flag);

	      for (j = 0; j < sa[0] * sa[1]; j++)
		{
		  cinf[j] = einf[j];
		  csup[j] = esup[j];
		}
	      free (dinf);
	      free (dsup);
	      free (einf);
	      free (esup);
	      free (finf);
	      free (fsup);
	      free (ginf);
	      free (gsup);
	      return 0;
	    }

	  // initialisation of f to the square or not
	  if (Bin[1] == 1)
	    {
	      for (j = 0; j < sa[0] * sa[1]; j++)
		{
		  finf[j] = dinf[j];
		  fsup[j] = dsup[j];
		}
	      Init_check = 1;
	    }

	  else
	    {
	      Init_check = 0;
	    }

	  for (i = 2; i < K; i++)
	    {
	      // Operation with squares
	      biasmulII4sci (dinf, dsup, sa, dinf, dsup, sa, einf, esup,
			     flag);

	      if (Bin[i] == 1)
		{
		  if (Init_check == 0)
		    {
		      for (j = 0; j < sa[0] * sa[1]; j++)
			{
			  finf[j] = einf[j];
			  fsup[j] = esup[j];
			}
		      Init_check = 1;
		    }

		  else
		    {
		      biasmulII4sci (finf, fsup, sa, einf, esup, sa, ginf,
				     gsup, flag);

		      for (j = 0; j < sa[0] * sa[1]; j++)
			{
			  finf[j] = ginf[j];
			  fsup[j] = gsup[j];
			}
		    }
		}

	      // Passing d to next square power
	      for (j = 0; j < sa[0] * sa[1]; j++)
		{
		  dinf[j] = einf[j];
		  dsup[j] = esup[j];
		}
	    }

	  // Finalisation in function of last binary value
	  if (Bin[0] == 1)
	    {
	      biasmulII4sci (finf, fsup, sa, ainf, asup, sa, cinf, csup,
			     flag);
	      free (dinf);
	      free (dsup);
	      free (einf);
	      free (esup);
	      free (finf);
	      free (fsup);
	      free (ginf);
	      free (gsup);
	      return 0;
	    }

	  else
	    {
	      for (i = 0; i < sa[0] * sa[1]; i++)
		{
		  cinf[i] = finf[i];
		  csup[i] = fsup[i];
		}
	      free (dinf);
	      free (dsup);
	      free (einf);
	      free (esup);
	      free (finf);
	      free (fsup);
	      free (ginf);
	      free (gsup);
	      return 0;
	    }
	}
    }
}
int
biasRpowerNI4sci (const double *ainf, const double *asup, const int *sa,
		  const int *n, double *cinf, double *csup, int *flag)
{
  int i;

  BIASINTERVAL *pA;
  BIASINTERVAL *pC;

  BiasFuncInit ();
  pA = (BIASINTERVAL *) malloc (sa[0] * sa[1] * sizeof (BIASINTERVAL));
  pC = (BIASINTERVAL *) malloc (sa[0] * sa[1] * sizeof (BIASINTERVAL));
  See (pA, pC);
  See (pC, pA);
  if (n[0] == 2)
    { 
      for (i = 0; i < sa[0] * sa[1]; i++)
	{
	  pA[i].inf = ainf[i];
	  pA[i].sup = asup[i];
	  BiasSqr (&pC[i], &(pA[i]));
	  cinf[i] = pC[i].inf;
	  csup[i] = pC[i].sup;
	}
      *flag = 0;
      free (pA);
      free (pC);
      return 0;
    }

  else
    { 
      for (i = 0; i < sa[0] * sa[1]; i++)
	{
	  pA[i].inf = ainf[i];
	  pA[i].sup = asup[i];
	  BiasPowerN (&pC[i], &(pA[i]), n[0]);
	  cinf[i] = pC[i].inf;
	  csup[i] = pC[i].sup;
	}
       *flag = 0;
       free (pA);
       free (pC);
       return 0;
    }
}



int
biasabsI4sci (const double *ainf, const double *asup, const int *sa,
	      double *cinf, double *csup, int *flag)
{
  //feclearexcept (FE_ALL_EXCEPT);
  int i;

  BIASINTERVAL *pA;
  BIASINTERVAL *pC;

  BiasFuncInit ();
  pA = (BIASINTERVAL *) malloc (sa[0] * sa[1] * sizeof (BIASINTERVAL));
  pC = (BIASINTERVAL *) malloc (sa[0] * sa[1] * sizeof (BIASINTERVAL));
  See (pA, pC);
  See (pC, pA);

  for (i = 0; i < sa[0] * sa[1]; i++)
    {
      pA[i].inf = ainf[i];
      pA[i].sup = asup[i];
      BiasIAbs (&pC[i], &(pA[i]));
      cinf[i] = pC[i].inf;
      csup[i] = pC[i].sup;
    }
  free (pA);
  free (pC);
  *flag = 0;
  return 0;
}

int
biasexpI4sci (const double *ainf, const double *asup, const int *sa,
	      double *cinf, double *csup, int *flag)
{
  BIASINTERVAL *pA;
  BIASINTERVAL *pC;

  int i;

  BiasFuncInit ();

  pA = (BIASINTERVAL *) malloc (sa[0] * sa[1] * sizeof (BIASINTERVAL));
  pC = (BIASINTERVAL *) malloc (sa[0] * sa[1] * sizeof (BIASINTERVAL));
  See (pA, pC);
  See (pC, pA);

  for (i = 0; i < sa[0] * sa[1]; i++)
    {
      pA[i].inf = ainf[i];
      pA[i].sup = asup[i];
      BiasExp (&pC[i], &(pA[i]));
      cinf[i] = pC[i].inf;
      csup[i] = pC[i].sup;
    }

  free (pA);
  free (pC);
  *flag = 0;
  return 0;
}

int
biaslogI4sci (const double *ainf, const double *asup, const int *sa,
	      double *cinf, double *csup, int *flag)
{
  double zero = 0.0;
  int i;

  BIASINTERVAL *pA;
  BIASINTERVAL *pC;

  BiasFuncInit ();

  pA = (BIASINTERVAL *) malloc (sa[0] * sa[1] * sizeof (BIASINTERVAL));
  pC = (BIASINTERVAL *) malloc (sa[0] * sa[1] * sizeof (BIASINTERVAL));
  See (pA, pC);
  See (pC, pA);

  for (i = 0; i < sa[0] * sa[1]; i++)
    {
      if (asup[i] < 0)
	{
	  *flag = -1;
	  return 0;
	}

      pA[i].inf = ainf[i];
      pA[i].sup = asup[i];

      if (BiasInR (&zero, &(pA[i])))
	{
	  cinf[i] = -HUGE_VAL;
	  csup[i] = log (pA[i].sup);
	}

      else
	{
	  BiasLog (&pC[i], &(pA[i]));
	  cinf[i] = pC[i].inf;
	  csup[i] = pC[i].sup;
	}
    }

  free (pA);
  free (pC);
  *flag = 0;
  return 0;
}

int
biaslog10I4sci (const double *ainf, const double *asup, const int *sa,
		double *cinf, double *csup, int *flag)
{
  int i;
  double zero = 0.0;

  BIASINTERVAL *pA;
  BIASINTERVAL *pC;

  BiasFuncInit ();
  pA = (BIASINTERVAL *) malloc (sa[0] * sa[1] * sizeof (BIASINTERVAL));
  pC = (BIASINTERVAL *) malloc (sa[0] * sa[1] * sizeof (BIASINTERVAL));
  See (pA, pC);
  See (pC, pA);

  for (i = 0; i < sa[0] * sa[1]; i++)
    {
      if (asup[i] < 0)
	{
	  *flag = -1;
	  return 0;
	}

      pA[i].inf = ainf[i];
      pA[i].sup = asup[i];

      if (BiasInR (&zero, &(pA[i])))
	{
	  cinf[i] = -HUGE_VAL;
	  csup[i] = log10 (pA[i].sup);
	}

      else
	{
	  BiasLog10 (&pC[i], &(pA[i]));
	  cinf[i] = pC[i].inf;
	  csup[i] = pC[i].sup;
	}
    }

  free (pA);
  free (pC);
  *flag = 0;
  return 0;
}

int
biassinI4sci (const double *ainf, const double *asup, const int *sa,
	      double *cinf, double *csup, int *flag)
{
  BIASINTERVAL *pA;
  BIASINTERVAL *pC;
  int i;

  BiasFuncInit ();
  pA = (BIASINTERVAL *) malloc (sa[0] * sa[1] * sizeof (BIASINTERVAL));
  pC = (BIASINTERVAL *) malloc (sa[0] * sa[1] * sizeof (BIASINTERVAL));
  See (pA, pC);
  See (pC, pA);

  for (i = 0; i < sa[0] * sa[1]; i++)
    {
      pA[i].inf = ainf[i];
      pA[i].sup = asup[i];
      BiasSin (&(pC[i]), &(pA[i]));
      cinf[i] = pC[i].inf;
      csup[i] = pC[i].sup;
    }

  free (pA);
  free (pC);
  *flag = 0;
  return 0;
}

int
biascosI4sci (const double *ainf, const double *asup, const int *sa,
	      double *cinf, double *csup, int *flag)
{
  BIASINTERVAL *pA;
  BIASINTERVAL *pC;
  int i;

  BiasFuncInit ();

  pA = (BIASINTERVAL *) malloc (sa[0] * sa[1] * sizeof (BIASINTERVAL));
  pC = (BIASINTERVAL *) malloc (sa[0] * sa[1] * sizeof (BIASINTERVAL));
  See (pA, pC);
  See (pC, pA);

  for (i = 0; i < sa[0] * sa[1]; i++)
    {
      pA[i].inf = ainf[i];
      pA[i].sup = asup[i];
      BiasCos (&pC[i], &(pA[i]));
      cinf[i] = pC[i].inf;
      csup[i] = pC[i].sup;
    }

  free (pA);
  free (pC);
  *flag = 0;
  return 0;
}

int
biastanI4sci (const double *ainf, const double *asup, const int *sa,
	      double *cinf, double *csup, int *flag)
{
  BIASINTERVAL *pA;
  BIASINTERVAL *pC;
  int i;
  BiasFuncInit ();
  pA = (BIASINTERVAL *) malloc (sa[0] * sa[1] * sizeof (BIASINTERVAL));
  pC = (BIASINTERVAL *) malloc (sa[0] * sa[1] * sizeof (BIASINTERVAL));
  See (pA, pC);
  See (pC, pA);

  for (i = 0; i < sa[0] * sa[1]; i++)
    {
      pA[i].inf = ainf[i];
      pA[i].sup = asup[i];
      BiasTan (&pC[i], &(pA[i]));
      cinf[i] = pC[i].inf;
      csup[i] = pC[i].sup;
    }
  free (pA);
  free (pC);
  *flag = 0;
  return 0;
}

int
biascotI4sci (const double *ainf, const double *asup, const int *sa,
	      double *cinf, double *csup, int *flag)
{
  BIASINTERVAL *pA;
  BIASINTERVAL *pC;
  int i;

  BiasFuncInit ();

  pA = (BIASINTERVAL *) malloc (sa[0] * sa[1] * sizeof (BIASINTERVAL));
  pC = (BIASINTERVAL *) malloc (sa[0] * sa[1] * sizeof (BIASINTERVAL));
  See (pA, pC);
  See (pC, pA);

  for (i = 0; i < sa[0] * sa[1]; i++)
    {
      pA[i].inf = ainf[i];
      pA[i].sup = asup[i];
      BiasCot (&pC[i], &(pA[i]));
      cinf[i] = pC[i].inf;
      csup[i] = pC[i].sup;
    }

  free (pA);
  free (pC);
  *flag = 0;
  return 0;
}


int
biasarcsinI4sci (const double *ainf, const double *asup, const int *sa,
		 double *cinf, double *csup, int *flag)
{
  BIASINTERVAL *pA;
  BIASINTERVAL *pC;
  int i;

  BiasFuncInit ();
  pA = (BIASINTERVAL *) malloc (sa[0] * sa[1] * sizeof (BIASINTERVAL));
  pC = (BIASINTERVAL *) malloc (sa[0] * sa[1] * sizeof (BIASINTERVAL));
  See (pA, pC);
  See (pC, pA);

  for (i = 0; i < sa[0] * sa[1]; i++)
    {
      pA[i].inf = ainf[i];
      pA[i].sup = asup[i];
      BiasArcSin (&pC[i], &(pA[i]));
      cinf[i] = pC[i].inf;
      csup[i] = pC[i].sup;
    }
  free (pA);
  free (pC);
  *flag = 0;
  return 0;
}

int
biasarccosI4sci (const double *ainf, const double *asup, const int *sa,
		 double *cinf, double *csup, int *flag)
{
  BIASINTERVAL *pA;
  BIASINTERVAL *pC;
  int i;

  BiasFuncInit ();

  pA = (BIASINTERVAL *) malloc (sa[0] * sa[1] * sizeof (BIASINTERVAL));
  pC = (BIASINTERVAL *) malloc (sa[0] * sa[1] * sizeof (BIASINTERVAL));
  See (pA, pC);
  See (pC, pA);

  for (i = 0; i < sa[0] * sa[1]; i++)
    {
      pA[i].inf = ainf[i];
      pA[i].sup = asup[i];
      BiasArcCos (&pC[i], &(pA[i]));
      cinf[i] = pC[i].inf;
      csup[i] = pC[i].sup;
    }

  free (pA);
  free (pC);
  *flag = 0;
  return 0;
}

int
biasarctanI4sci (const double *ainf, const double *asup, const int *sa,
		 double *cinf, double *csup, int *flag)
{
  BIASINTERVAL *pA;
  BIASINTERVAL *pC;
  int i;

  BiasFuncInit ();

  pA = (BIASINTERVAL *) malloc (sa[0] * sa[1] * sizeof (BIASINTERVAL));
  pC = (BIASINTERVAL *) malloc (sa[0] * sa[1] * sizeof (BIASINTERVAL));
  See (pA, pC);
  See (pC, pA);

  for (i = 0; i < sa[0] * sa[1]; i++)
    {
      pA[i].inf = ainf[i];
      pA[i].sup = asup[i];
      BiasArcTan (&pC[i], &(pA[i]));
      cinf[i] = pC[i].inf;
      csup[i] = pC[i].sup;
    }

  free (pA);
  free (pC);
  *flag = 0;
  return 0;
}

int
biasarccotI4sci (const double *ainf, const double *asup, const int *sa,
		 double *cinf, double *csup, int *flag)
{
  BIASINTERVAL *pA;
  BIASINTERVAL *pC;
  int i;

  BiasFuncInit ();

  pA = (BIASINTERVAL *) malloc (sa[0] * sa[1] * sizeof (BIASINTERVAL));
  pC = (BIASINTERVAL *) malloc (sa[0] * sa[1] * sizeof (BIASINTERVAL));
  See (pA, pC);
  See (pC, pA);

  for (i = 0; i < sa[0] * sa[1]; i++)
    {
      pA[i].inf = ainf[i];
      pA[i].sup = asup[i];
      BiasArcCot (&pC[i], &(pA[i]));
      cinf[i] = pC[i].inf;
      csup[i] = pC[i].sup;
    }

  free (pA);
  free (pC);
  *flag = 0;
  return 0;
}

int
biassinhI4sci (const double *ainf, const double *asup, const int *sa,
	       double *cinf, double *csup, int *flag)
{
  BIASINTERVAL *pA;
  BIASINTERVAL *pC;
  int i;

  BiasFuncInit ();

  pA = (BIASINTERVAL *) malloc (sa[0] * sa[1] * sizeof (BIASINTERVAL));
  pC = (BIASINTERVAL *) malloc (sa[0] * sa[1] * sizeof (BIASINTERVAL));
  See (pA, pC);
  See (pC, pA);

  for (i = 0; i < sa[0] * sa[1]; i++)
    {
      pA[i].inf = ainf[i];
      pA[i].sup = asup[i];
      BiasSinh (&pC[i], &(pA[i]));
      cinf[i] = pC[i].inf;
      csup[i] = pC[i].sup;
    }

  free (pA);
  free (pC);
  *flag = 0;
  return 0;
}

int
biascoshI4sci (const double *ainf, const double *asup, const int *sa,
	       double *cinf, double *csup, int *flag)
{
  BIASINTERVAL *pA;
  BIASINTERVAL *pC;
  int i;

  BiasFuncInit ();
  pA = (BIASINTERVAL *) malloc (sa[0] * sa[1] * sizeof (BIASINTERVAL));
  pC = (BIASINTERVAL *) malloc (sa[0] * sa[1] * sizeof (BIASINTERVAL));
  See (pA, pC);
  See (pC, pA);

  for (i = 0; i < sa[0] * sa[1]; i++)
    {
      pA[i].inf = ainf[i];
      pA[i].sup = asup[i];
      BiasCosh (&pC[i], &(pA[i]));
      cinf[i] = pC[i].inf;
      csup[i] = pC[i].sup;
    }

  free (pA);
  free (pC);
  *flag = 0;
  return 0;
}

int
biastanhI4sci (const double *ainf, const double *asup, const int *sa,
	       double *cinf, double *csup, int *flag)
{
  BIASINTERVAL *pA;
  BIASINTERVAL *pC;
  int i;

  BiasFuncInit ();
  pA = (BIASINTERVAL *) malloc (sa[0] * sa[1] * sizeof (BIASINTERVAL));
  pC = (BIASINTERVAL *) malloc (sa[0] * sa[1] * sizeof (BIASINTERVAL));
  See (pA, pC);
  See (pC, pA);

  for (i = 0; i < sa[0] * sa[1]; i++)
    {
      pA[i].inf = ainf[i];
      pA[i].sup = asup[i];
      BiasTanh (&pC[i], &(pA[i]));
      cinf[i] = pC[i].inf;
      csup[i] = pC[i].sup;
    }

  free (pA);
  free (pC);
  *flag = 0;
  return 0;
}

int
biascothI4sci (const double *ainf, const double *asup, const int *sa,
	       double *cinf, double *csup, int *flag)
{
  BIASINTERVAL *pA;
  BIASINTERVAL *pC;
  int i;

  BiasFuncInit ();

  pA = (BIASINTERVAL *) malloc (sa[0] * sa[1] * sizeof (BIASINTERVAL));
  pC = (BIASINTERVAL *) malloc (sa[0] * sa[1] * sizeof (BIASINTERVAL));
  See (pA, pC);
  See (pC, pA);

  for (i = 0; i < sa[0] * sa[1]; i++)
    {
      pA[i].inf = ainf[i];
      pA[i].sup = asup[i];
      BiasCoth (&pC[i], &(pA[i]));
      cinf[i] = pC[i].inf;
      csup[i] = pC[i].sup;
    }

  free (pA);
  free (pC);
  *flag = 0;
  return 0;
}

int
biasarsinhI4sci (const double *ainf, const double *asup, const int *sa,
		 double *cinf, double *csup, int *flag)
{
  BIASINTERVAL *pA;
  BIASINTERVAL *pC;
  int i;

  BiasFuncInit ();

  pA = (BIASINTERVAL *) malloc (sa[0] * sa[1] * sizeof (BIASINTERVAL));
  pC = (BIASINTERVAL *) malloc (sa[0] * sa[1] * sizeof (BIASINTERVAL));
  See (pA, pC);
  See (pC, pA);

  for (i = 0; i < sa[0] * sa[1]; i++)
    {
      pA[i].inf = ainf[i];
      pA[i].sup = asup[i];
      BiasArSinh (&pC[i], &(pA[i]));
      cinf[i] = pC[i].inf;
      csup[i] = pC[i].sup;
    }

  free (pA);
  free (pC);
  *flag = 0;
  return 0;
}

int
biasarcoshI4sci (const double *ainf, const double *asup, const int *sa,
		 double *cinf, double *csup, int *flag)
{
  BIASINTERVAL *pA;
  BIASINTERVAL *pC;
  int i;

  BiasFuncInit ();

  pA = (BIASINTERVAL *) malloc (sa[0] * sa[1] * sizeof (BIASINTERVAL));
  pC = (BIASINTERVAL *) malloc (sa[0] * sa[1] * sizeof (BIASINTERVAL));
  See (pA, pC);
  See (pC, pA);

  for (i = 0; i < sa[0] * sa[1]; i++)
    {
      pA[i].inf = ainf[i];
      pA[i].sup = asup[i];
      BiasArCosh (&pC[i], &(pA[i]));
      cinf[i] = pC[i].inf;
      csup[i] = pC[i].sup;
    }

  free (pA);
  free (pC);
  *flag = 0;
  return 0;
}

int
biasartanhI4sci (const double *ainf, const double *asup, const int *sa,
		 double *cinf, double *csup, int *flag)
{
  BIASINTERVAL *pA;
  BIASINTERVAL *pC;
  int i;

  BiasFuncInit ();
  pA = (BIASINTERVAL *) malloc (sa[0] * sa[1] * sizeof (BIASINTERVAL));
  pC = (BIASINTERVAL *) malloc (sa[0] * sa[1] * sizeof (BIASINTERVAL));
  See (pA, pC);
  See (pC, pA);

  for (i = 0; i < sa[0] * sa[1]; i++)
    {
      pA[i].inf = ainf[i];
      pA[i].sup = asup[i];
      BiasArTanh (&pC[i], &(pA[i]));
      cinf[i] = pC[i].inf;
      csup[i] = pC[i].sup;
    }

  free (pA);
  free (pC);
  *flag = 0;
  return 0;
}

int
biasarcothI4sci (const double *ainf, const double *asup, const int *sa,
		 double *cinf, double *csup, int *flag)
{
  BIASINTERVAL *pA;
  BIASINTERVAL *pC;
  int i;

  BiasFuncInit ();

  pA = (BIASINTERVAL *) malloc (sa[0] * sa[1] * sizeof (BIASINTERVAL));
  pC = (BIASINTERVAL *) malloc (sa[0] * sa[1] * sizeof (BIASINTERVAL));
  See (pA, pC);
  See (pC, pA);

  for (i = 0; i < sa[0] * sa[1]; i++)
    {
      pA[i].inf = ainf[i];
      pA[i].sup = asup[i];
      BiasArCoth (&pC[i], &(pA[i]));
      cinf[i] = pC[i].inf;
      csup[i] = pC[i].sup;
    }

  free (pA);
  free (pC);
  *flag = 0;
  return 0;
}

int
checkmul4sci (double *a, double *b, double *c, double *d, int *t)
{
  t[0] = 1;

  if ((a[0] == 0 && b[0] == 0) || (c[0] == 0 && d[0] == 0))
    {
      a[0] = 0;
      b[0] = 0;
      c[0] = 1;
      d[0] = 1;
      t[0] = 0;
      return 1;
    }

  if (a[0] == -HUGE_VAL && b[0] == HUGE_VAL)

    {
      if (c[0] == 0 && d[0] == 0)
	{
	  a[0] = 0;
	  b[0] = 0;
	  c[0] = 1;
	  d[0] = 1;
	  t[0] = 1;
	  return 1;
	}
      else
	{
	  a[0] = -HUGE_VAL;
	  b[0] = HUGE_VAL;
	  c[0] = 1;
	  d[0] = 2;
	  t[0] = 0;
	  return 1;
	}
    }

  if (a[0] == HUGE_VAL && b[0] == HUGE_VAL)

    {
      if (d[0] == 0)

	{
	  a[0] = -HUGE_VAL;
	  b[0] = 0;
	  c[0] = 1;
	  d[0] = 1;
	  t[0] = 0;
	  return 1;
	}

      else if (c[0] == 0)

	{
	  a[0] = 0;
	  b[0] = HUGE_VAL;
	  c[0] = 1;
	  d[0] = 1;
	  t[0] = 0;
	  return 1;
	}
    }

  else if (a[0] == -HUGE_VAL && b[0] == -HUGE_VAL)

    {
      if (d[0] == 0)

	{
	  a[0] = 0;
	  b[0] = HUGE_VAL;
	  c[0] = 1;
	  d[0] = 1;
	  t[0] = 0;
	  return 1;
	}

      else if (c[0] == 0)

	{
	  a[0] = -HUGE_VAL;
	  b[0] = 0;
	  c[0] = 1;
	  d[0] = 1;
	  t[0] = 0;
	  return 1;
	}
    }

  else if (a[0] == -HUGE_VAL && d[0] == 0)

    {
      if (b[0] > 0)

	{
	  a[0] = c[0] * b[0];
	  b[0] = HUGE_VAL;
	  c[0] = 1;
	  d[0] = 1;
	  t[0] = 0;
	  return 1;
	}

      else

	{
	  a[0] = 0;
	  b[0] = HUGE_VAL;
	  c[0] = 1;
	  d[0] = 1;
	  t[0] = 0;
	  return 1;
	}
    }

  else if (b[0] == HUGE_VAL && d[0] == 0)

    {
      if (a[0] < 0)

	{
	  b[0] = a[0] * c[0];
	  a[0] = -HUGE_VAL;
	  c[0] = 1;
	  d[0] = 1;
	  t[0] = 0;
	  return 1;
	}

      else

	{
	  b[0] = 0;
	  a[0] = -HUGE_VAL;
	  c[0] = 1;
	  d[0] = 1;
	  t[0] = 0;
	  return 1;
	}
    }
}

int
DIYMulMIMI (BIASINTERVAL * pC, BIASINTERVAL * pA, BIASINTERVAL * pB,
	    const int *sa, const int *sb)
{

  int i, j, k;
  int t[1];
  BIASINTERVAL *pF;
  BIASINTERVAL pX[1];
  BIASINTERVAL pY[1];
  double ai[1];
  double as[1];
  double bi[1];
  double bs[1];

  pF = (BIASINTERVAL *) malloc (sa[0] * sb[1] * sizeof (BIASINTERVAL));

  for (i = 0; i < sa[0]; i++)
    {
      for (k = 0; k < sb[1]; k++)
	{
	  pF[k + i * sb[1]].inf = 0;
	  pF[k + i * sb[1]].sup = 0;

	  for (j = 0; j < sa[1]; j++)
	    {
	      ai[0] = pA[j + i * sa[1]].inf;
	      as[0] = pA[j + i * sa[1]].sup;
	      bi[0] = pB[k + j * sb[1]].inf;
	      bs[0] = pB[k + j * sb[1]].sup;

	      checkmul4sci (ai, as, bi, bs, t);

	      if (t[0] != 0)
		{
		  checkmul4sci (bi, bs, ai, as, t);
		}

	      pX[0].inf = ai[0];
	      pX[0].sup = as[0];
	      pY[0].inf = bi[0];
	      pY[0].sup = bs[0];

	      BiasMacII (&(pF[k + i * sb[1]]), pX, pY);
	    }

	  pC[k + i * sb[1]] = pF[k + i * sb[1]];
	}
    }
  free (pF);
  return 1;
}


int
DIYMulIMI (BIASINTERVAL * pC, BIASINTERVAL * pA, BIASINTERVAL * pB,
	   const int *sb)
{
  int i, k;
  int t[1];
  BIASINTERVAL *pF;
  BIASINTERVAL pX, pY;
  double ai[1];
  double as[1];
  double bi[1];
  double bs[1];

  pF = (BIASINTERVAL *) malloc (sb[0] * sb[1] * sizeof (BIASINTERVAL));

  for (i = 0; i < sb[0]; i++)
    {
      for (k = 0; k < sb[1]; k++)
	{
	  pF[k + i * sb[1]].inf = 0;
	  pF[k + i * sb[1]].sup = 0;

	  *ai = pA[0].inf;
	  *as = pA[0].sup;
	  *bi = pB[k + i * sb[1]].inf;
	  *bs = pB[k + i * sb[1]].sup;

	  checkmul4sci (ai, as, bi, bs, t);

	  if (*t != 0)
	    {
	      checkmul4sci (bi, bs, ai, as, t);
	    }

	  pX.inf = *ai;
	  pX.sup = *as;
	  pY.inf = *bi;
	  pY.sup = *bs;

	  BiasMacII (&(pF[k + i * sb[1]]), &(pX), &(pY));


	  pC[k + i * sb[1]] = pF[k + i * sb[1]];
	}
    }

  free (pF);
  return 1;
}


int
DIYMulMRMI (BIASINTERVAL * pC, const double *A, BIASINTERVAL * pB,
	    const int *sa, const int *sb)
{
  int i, j, k;
  int t[1];
  BIASINTERVAL *pF;
  BIASINTERVAL pX, pY;
  double ai[1];
  double as[1];
  double bi[1];
  double bs[1];

  pF = (BIASINTERVAL *) malloc (sa[0] * sb[1] * sizeof (BIASINTERVAL));

  for (i = 0; i < sa[0]; i++)
    {
      for (k = 0; k < sb[1]; k++)
	{
	  pF[k + i * sb[1]].inf = 0;
	  pF[k + i * sb[1]].sup = 0;

	  for (j = 0; j < sa[1]; j++)
	    {
	      *ai = A[j + i * sa[1]];
	      *as = A[j + i * sa[1]];
	      *bi = pB[k + j * sb[1]].inf;
	      *bs = pB[k + j * sb[1]].sup;

	      checkmul4sci (ai, as, bi, bs, t);

	      if (*t != 0)
		{
		  checkmul4sci (bi, bs, ai, as, t);
		}

	      pX.inf = *ai;
	      pX.sup = *as;
	      pY.inf = *bi;
	      pY.sup = *bs;

	      BiasMacII (&(pF[k + i * sb[1]]), &(pX), &(pY));
	    }

	  pC[k + i * sb[1]] = pF[k + i * sb[1]];
	}
    }

  free (pF);
  return 1;
}

int
DIYMulMIMR (BIASINTERVAL * pC, BIASINTERVAL * pA, const double  * B,
	    const int *sa, const int *sb)
{
  int i, j, k;
  int t[1];
  BIASINTERVAL *pF;
  BIASINTERVAL pX, pY;
  double ai[1];
  double as[1];
  double bi[1];
  double bs[1];

  pF = (BIASINTERVAL *) malloc (sa[0] * sb[1] * sizeof (BIASINTERVAL));

  for (i = 0; i < sa[0]; i++)
    {
      for (k = 0; k < sb[1]; k++)
	{
	  pF[k + i * sb[1]].inf = 0;
	  pF[k + i * sb[1]].sup = 0;

	  for (j = 0; j < sa[1]; j++)
	    {
	      *ai = pA[j + i * sa[1]].inf;
	      *as = pA[j + i * sa[1]].sup;
	      *bi = B[k + j * sb[1]];
	      *bs = B[k + j * sb[1]];

	      checkmul4sci (ai, as, bi, bs, t);

	      if (*t != 0)
		{
		  checkmul4sci (bi, bs, ai, as, t);
		}

	      pX.inf = *ai;
	      pX.sup = *as;
	      pY.inf = *bi;
	      pY.sup = *bs;

	      BiasMacII (&(pF[k + i * sb[1]]), &(pX), &(pY));
	    }

	  pC[k + i * sb[1]] = pF[k + i * sb[1]];
	}
    }

  free (pF);
  return 1;
}

int
DIYMulIMR (BIASINTERVAL * pC, BIASINTERVAL * pB, const double *A,
	   const int *sa)
{
  int i, k;
  int t[1];
  BIASINTERVAL *pF;
  BIASINTERVAL pX, pY;
  double ai[1];
  double as[1];
  double bi[1];
  double bs[1];


  pF = (BIASINTERVAL *) malloc (sa[0] * sa[1] * sizeof (BIASINTERVAL));

  for (i = 0; i < sa[0]; i++)
    {
      for (k = 0; k < sa[1]; k++)
	{
	  pF[k + i * sa[1]].inf = 0;
	  pF[k + i * sa[1]].sup = 0;

	  *ai = A[k + i * sa[1]];
	  *as = A[k + i * sa[1]];
	  *bi = pB[0].inf;
	  *bs = pB[0].sup;

	  checkmul4sci (ai, as, bi, bs, t);

	  if (*t != 0)
	    {
	      checkmul4sci (bi, bs, ai, as, t);
	    }

	  pX.inf = *ai;
	  pX.sup = *as;
	  pY.inf = *bi;
	  pY.sup = *bs;

	  BiasMacII (&(pF[k + i * sa[1]]), &(pX), &(pY));

	  pC[k + i * sa[1]] = pF[k + i * sa[1]];
	}
    }

  free (pF);
  return 1;
}

int
DIYMulRMI (BIASINTERVAL * pC, const double *A, BIASINTERVAL * pB,
	   const int *sb)
{
  int i, k;
  int t[1];
  BIASINTERVAL *pF;
  BIASINTERVAL pX, pY;
  double ai[1];
  double as[1];
  double bi[1];
  double bs[1];


  pF = (BIASINTERVAL *) malloc (sb[0] * sb[1] * sizeof (BIASINTERVAL));

  for (i = 0; i < sb[0]; i++)
    {
      for (k = 0; k < sb[1]; k++)
	{
	  pF[k + i * sb[1]].inf = 0;
	  pF[k + i * sb[1]].sup = 0;

	  *ai = A[0];
	  *as = A[0];
	  *bi = pB[k + i * sb[1]].inf;
	  *bs = pB[k + i * sb[1]].sup;

	  checkmul4sci (ai, as, bi, bs, t);

	  if (*t != 0)
	    {
	      checkmul4sci (bi, bs, ai, as, t);
	    }

	  pX.inf = *ai;
	  pX.sup = *as;
	  pY.inf = *bi;
	  pY.sup = *bs;

	  BiasMacII (&(pF[k + i * sb[1]]), &(pX), &(pY));

	  pC[k + i * sb[1]] = pF[k + i * sb[1]];
	}
    }

  free (pF);
  return 1;
}

int
DIYKMulMIMI (BIASINTERVAL * pC, BIASINTERVAL * pA, BIASINTERVAL * pB,
	     const int *sa, const int *sb)
{
  int i, j;
  BIASINTERVAL *pE;

  for (i = 0; i < sa[0] * sa[1]; i++)
    {
      pE = (BIASINTERVAL *) malloc (sb[0] * sb[1] * sizeof (BIASINTERVAL));

      DIYMulIMI (pE, &(pA[i]), pB, sb);

      for (j = 0; j < sb[0] * sb[1]; j++)
	{
	  pC[j % sb[1] + (i % sa[1]) * sb[1]
	     + (j / sb[1] + (i / sa[1]) * sb[0]) * sa[1] * sb[1]] = pE[j];

	}
      free (pE);
    }
  return 1;
}

int
DIYKMulMIMR (BIASINTERVAL * pC, BIASINTERVAL * pA, double *rB, const int *sa,
	     const int *sb)
{
  int i, j;
  BIASINTERVAL *pE;

  for (i = 0; i < sa[0] * sa[1]; i++)
    {
      pE = (BIASINTERVAL *) malloc (sb[0] * sb[1] * sizeof (BIASINTERVAL));

      DIYMulIMR (pE, &(pA[i]), rB, sb);

      for (j = 0; j < sb[0] * sb[1]; j++)
	{
	  pC[j % sb[1] + (i % sa[1]) * sb[1]
	     + (j / sb[1] + (i / sa[1]) * sb[0]) * sa[1] * sb[1]] = pE[j];
	}
      free (pE);
    }
  return 1;
}

int
DIYKMulMRMI (BIASINTERVAL * pC, double *rA, BIASINTERVAL * pB, const int *sa,
	     const int *sb)
{
  int i, j;
  BIASINTERVAL *pE;

  for (i = 0; i < sa[0] * sa[1]; i++)
    {
      pE = (BIASINTERVAL *) malloc (sb[0] * sb[1] * sizeof (BIASINTERVAL));

      DIYMulRMI (pE, &(rA[i]), pB, sb);

      for (j = 0; j < sb[0] * sb[1]; j++)
	{
	  pC[j % sb[1] + (i % sa[1]) * sb[1]
	     + (j / sb[1] + (i / sa[1]) * sb[0]) * sa[1] * sb[1]] = pE[j];

	}
      free (pE);
    }
  return 1;
}

int
DIYDivII (BIASINTERVAL * pC, BIASINTERVAL * pA, BIASINTERVAL * pB)
{
  if (pB->inf == 0)
    {
      if (pB->sup == HUGE_VAL)
	{
	  if (pA->sup <= 0)
	    {
	      pC->inf = -HUGE_VAL;
	      pC->sup = 0;
	      return 1;
	    }
	  else if (pA->inf >= 0)
	    {
	      pC->inf = 0;
	      pC->sup = HUGE_VAL;
	      return 1;
	    }
	  else
	    {
	      pC->inf = -HUGE_VAL;
	      pC->sup = HUGE_VAL;
	      return 1;
	    }
	}
      else
	{
	  if (pA->sup <= 0)
	    {
	      pC->inf = -HUGE_VAL;
	      pC->sup = (pA->sup) / (pB->sup);
	      return 1;
	    }
	  else
	    {
	      pC->inf = (pA->inf) / (pB->sup);
	      pC->sup = HUGE_VAL;
	      return 1;
	    }
	}
    }

  if (fabs (pB->inf) == HUGE_VAL)
    {
      if (fabs (pB->sup) == HUGE_VAL)
	{
	  pC->inf = 0;
	  pC->sup = 0;
	  return 1;
	}
      else if ((pB->sup) < 0)
	{
	  if (pA->sup <= 0)
	    {
	      pC->inf = 0;
	      pC->sup = (pA->inf) / (pB->sup);
	      return 1;
	    }
	  else if (pA->inf >= 0)
	    {
	      pC->sup = 0;
	      pC->inf = (pA->sup) / (pB->sup);
	      return 1;
	    }
	  else
	    {
	      pC->sup = (pA->inf) / (pB->sup);
	      pC->inf = (pA->sup) / (pB->sup);
	      return 1;
	    }
	}
    }

  if (pB->sup == 0)
    {
      if (pB->inf == -HUGE_VAL)
	{
	  if (pA->sup <= 0)
	    {
	      pC->inf = 0;
	      pC->sup = HUGE_VAL;
	      return 1;
	    }
	  else if (pA->inf >= 0)
	    {
	      pC->inf = HUGE_VAL;
	      pC->sup = 0;
	      return 1;
	    }
	  else
	    {
	      pC->inf = -HUGE_VAL;
	      pC->sup = HUGE_VAL;
	      return 1;
	    }
	}
      else
	{
	  if (pA->sup <= 0)
	    {
	      pC->sup = HUGE_VAL;
	      pC->inf = (pA->sup) / (pB->inf);
	      return 1;
	    }
	  else
	    {
	      pC->inf = -HUGE_VAL;
	      pC->sup = (pA->inf) / (pB->inf);
	      return 1;
	    }
	}
    }

  if (fabs (pB->sup) == HUGE_VAL)
    {
      if (fabs (pB->inf) == HUGE_VAL)
	{
	  pC->inf = 0;
	  pC->sup = 0;
	  return 1;
	}
      else
	{
	  if (pA->sup <= 0)
	    {
	      pC->sup = 0;
	      pC->inf = (pA->inf) / (pB->inf);
	      return 1;
	    }
	  else if (pA->inf >= 0)
	    {
	      pC->inf = 0;
	      pC->sup = (pA->sup) / (pB->inf);
	      return 1;
	    }
	  else
	    {
	      pC->inf = (pA->inf) / (pB->inf);
	      pC->sup = (pA->sup) / (pB->inf);
	      return 1;
	    }
	}
    }
}

int
DIYDivRI (BIASINTERVAL * pC, double *rA, BIASINTERVAL * pB)
{
  double a = *rA;

  if (a == 0)
    {
      pC->inf = 0;
      pC->sup = 0;
      return 1;
    }
  else
    {
      BIASINTERVAL pA;
      pA.inf = a;
      pA.sup = a;
      DIYDivII (pC, &(pA), pB);
      return 1;
    }
}
