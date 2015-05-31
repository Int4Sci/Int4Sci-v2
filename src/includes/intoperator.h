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

#ifndef __OPERATE_H__
#define __OPERATE_H__

#if defined (__cplusplus)
extern "C"
{
	#endif

	int operate (const double *, const double *, const int *,
							 const double *, const double *, const int *,
							 char, double *, double *, int *);

	int power (const double *, const double *, const int *,
						 char, const int *, double *, double *, int *);

	int unary (const double *, const double *, const int *,
						 char, double *, double *, int *);

	int roundmode (const int *);
#if defined (__cplusplus)
}
#endif

#endif
