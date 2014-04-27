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


#ifndef __PROFIL_ILSR__H__
#define __PROFIL_ILSR__H__

#include "Utilities.h"
#include "Constants.h"
#include "Functions.h"
#include <fstream>
#include <stdio.h>
#include "IntervalVector.h"
#include "IntervalMatrix.h"
#include "LSS.h"

#include "IntegerVector.h"

/*------------------------------------------------------
 * Generic Declaration ils (interval linear solving) function
 *------------------------------------------------------*/
typedef int (*fctILS)(const INTERVAL_MATRIX& A, 
 		      const INTERVAL_VECTOR& b, 
 		      INTERVAL_VECTOR & x); 

// Preconditionning of ils functions
int pils(fctILS prog,const INTERVAL_MATRIX& A,const INTERVAL_VECTOR& b, INTERVAL_VECTOR &x);

// Contraction algorithm
// Krawczyck
int ilsk(const INTERVAL_MATRIX& A, const INTERVAL_VECTOR& b, INTERVAL_VECTOR &x);
int ilskgs(const INTERVAL_MATRIX& A, const INTERVAL_VECTOR& b, INTERVAL_VECTOR &x);
// Gauss elimination
int ilsge(const INTERVAL_MATRIX& A, const INTERVAL_VECTOR& b, INTERVAL_VECTOR & x);
// Bareiss
int ilsb(const INTERVAL_MATRIX& A, const INTERVAL_VECTOR& b, INTERVAL_VECTOR & x);
// Gauss sidel
int ilsgs(const INTERVAL_MATRIX& A, const INTERVAL_VECTOR& b, INTERVAL_VECTOR & x);
// Hansen Bliek
int  ilshb(const INTERVAL_MATRIX& A, const INTERVAL_VECTOR& b, INTERVAL_VECTOR & x);

REAL INorm(INTERVAL_VECTOR);
REAL INorm(INTERVAL_MATRIX);

int MatInverse(MATRIX, MATRIX &);


REAL Perimeter(INTERVAL_VECTOR); 

#endif
