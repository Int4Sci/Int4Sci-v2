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

// =================================//
// This function allows the resolution of F(x) = 0//
// for x interval with the newton method //

function y=F(x)
       y=(x-1)*(x-2)*(x-3);
endfunction

function y=Fp(x)
       y=(x-1)*(x-3)+(x-2)*(x-3)+(x-2)*(x-1);
endfunction

// Complete Newton

Zsol = [];
Zlst = interval(0,10);

accu=0;

while length(Zlst) > 0
  accu=accu+1;
  Ztmp = Zlst(1);
  Zlst= Zlst(2:$);
  div = extdiv(F(mid(Ztmp)),Fp(Ztmp));
  for i =1:length(div)
     Zint = intersection(Ztmp,-div(i)+mid(Ztmp));
     if (~isempty(Zint))
       if rad(Zint) < 1e-3
      Zsol = [Zsol Zint];
    else
      Zlst = [Zlst Zint];
    end
      end;
   end;
end
accu
Zsol