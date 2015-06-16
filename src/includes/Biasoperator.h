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


#ifndef __BIASOPERATOR_H__
#define __BIASOPERATOR_H__

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fenv.h>

#include <BIAS/Bias2.h>
#include <BIAS/BiasF.h>

#include <intoperator.h>
//#include <operate.h>


#ifndef max
#define max(X,Y) ((X)<(Y)? (Y):(X))
#endif

#define iterf2c(k,nl,nc) ((k%(nc))*(nl)+(k/(nc)))
#define iterc2f(k,nl,nc) ((k%(nl))*(nc)+(k/(nl)))

#if defined (__cplusplus)
extern "C"
{
#endif

  int biasaddII4sci (const double *, const double *, const int *,
		     const double *, const double *, const int *,
		     double *, double *, int *);

  int biasaddIR4sci (const double *, const double *, const int *,
		     const double *, const int *, double *, double *, int *);

  int biassubII4sci (const double *, const double *, const int *,
		     const double *, const double *, const int *,
		     double *, double *, int *);

  int biassubIR4sci (const double *, const double *, const int *,
		     const double *, const int *, double *, double *, int *);

  int biasmulII4sci (const double *, const double *, const int *,
		     const double *, const double *, const int *,
		     double *, double *, int *);

  int biasmulIR4sci (const double *, const double *, const int *,
		     const double *, const int *, double *, double *, int *);

  int biasmulRI4sci (const double *, const int *,
		     const double *, const double *, const int *,
		     double *, double *, int *);

  int biasRmulII4sci (const double *, const double *, const int *,
		      const double *, const double *, const int *,
		      double *, double *, int *);

  int biasRmulIR4sci (const double *, const double *, const int *,
		      const double *, const int *, double *, double *, int *);

  int biasRLmulII4sci (const double *, const double *, const int *,
		       const double *, const double *, const int *,
		       double *, double *, int *);

  int biasRLmulIR4sci (const double *, const double *, const int *,
		       const double *, const int *,
		       double *, double *, int *);

  int biasRLmulRI4sci (const double *, const int *,
		       const double *, const double *, const int *,
		       double *, double *, int *);

  int biasdivII4sci (const double *, const double *, const int *,
		     const double *, const double *, const int *,
		     double *, double *, int *);

  int biasdivRI4sci (const double *, const int *,
		     const double *, const double *, const int *,
		     double *, double *, int *);

  int biasdivIR4sci (const double *, const double *, const int *,
		     const double *, const int *, double *, double *, int *);

  int biasRdivII4sci (const double *, const double *, const int *,
		      const double *, const double *, const int *,
		      double *, double *, int *);

  int biasRdivRI4sci (const double *, const int *,
		      const double *, const double *, const int *,
		      double *, double *, int *);

  int biasRdivIR4sci (const double *, const double *, const int *,
		      const double *, const int *, double *, double *, int *);

  int biasrootNI4sci (const double *, const double *, const int *,
		      const int *, double *, double *, int *);

  int biasRpowerNI4sci (const double *, const double *, const int *,
			const int *, double *, double *, int *);

  int biaspowerNI4sci (const double *, const double *, const int *,
		       const int *, double *, double *, int *);

  int biasabsI4sci (const double *, const double *, const int *,
		    double *, double *, int *);

  int biasexpI4sci (const double *, const double *, const int *,
		    double *, double *, int *);

  int biaslogI4sci (const double *, const double *, const int *,
		    double *, double *, int *);

  int biaslog10I4sci (const double *, const double *, const int *,
		      double *, double *, int *);

  int biassinI4sci (const double *, const double *, const int *,
		    double *, double *, int *);

  int biascosI4sci (const double *, const double *, const int *,
		    double *, double *, int *);

  int biastanI4sci (const double *, const double *, const int *,
		    double *, double *, int *);

  int biascotI4sci (const double *, const double *, const int *,
		    double *, double *, int *);

  int biasarcsinI4sci (const double *, const double *, const int *,
		       double *, double *, int *);

  int biasarccosI4sci (const double *, const double *, const int *,
		       double *, double *, int *);

  int biasarctanI4sci (const double *, const double *, const int *,
		       double *, double *, int *);

  int biasarccotI4sci (const double *, const double *, const int *,
		       double *, double *, int *);

  int biassinhI4sci (const double *, const double *, const int *,
		     double *, double *, int *);

  int biascoshI4sci (const double *, const double *, const int *,
		     double *, double *, int *);

  int biastanhI4sci (const double *, const double *, const int *,
		     double *, double *, int *);

  int biascothI4sci (const double *, const double *, const int *,
		     double *, double *, int *);

  int biasarsinhI4sci (const double *, const double *, const int *,
		       double *, double *, int *);

  int biasarcoshI4sci (const double *, const double *, const int *,
		       double *, double *, int *);

  int biasartanhI4sci (const double *, const double *, const int *,
		       double *, double *, int *);

  int biasarcothI4sci (const double *, const double *, const int *,
		       double *, double *, int *);

  int DIYMulMIMI (BIASINTERVAL *, BIASINTERVAL *, BIASINTERVAL *,
		  const int[2], const int[2]);

  int DIYMulIMI (BIASINTERVAL *, BIASINTERVAL *, BIASINTERVAL *, const int *);

  int DIYMulMRMI (BIASINTERVAL *, const double *, BIASINTERVAL *, const int *,
		  const int *);

  int DIYMulMIMR (BIASINTERVAL *, BIASINTERVAL *, const double *, const int *,
		  const int *);

  int DIYMulRMI (BIASINTERVAL *, const double *, BIASINTERVAL *, const int *);

  int DIYMulIMR (BIASINTERVAL *, BIASINTERVAL *, const double *, const int *);

  int DIYKMulMIMI (BIASINTERVAL *, BIASINTERVAL *, BIASINTERVAL *,
		   const int *, const int *);

  int DIYKMulMIMR (BIASINTERVAL *, BIASINTERVAL *, double *, const int *,
		   const int *);

  int DIYKMulMRMI (BIASINTERVAL *, double *, BIASINTERVAL *, const int *,
		   const int *);

  int DIYDivII (BIASINTERVAL *, BIASINTERVAL *, BIASINTERVAL *);

  int DIYDivRI (BIASINTERVAL *, double *, BIASINTERVAL *);

  int checkmul4sci (double *, double *, double *, double *, int *);

#if defined (__cplusplus)
}
#endif

#endif
