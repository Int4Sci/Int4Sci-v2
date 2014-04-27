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

timer();

ee = "ERROR";

t = ee;
a = ee;
r = ee;
c = ee;
b = ee;
s = ee;


ierr=[];
disp('test for interval matrices insertion');
t = ones(10,10)

disp('insertion of intervals in t');
ierr(1) = execstr("t([1:2],[6:8]) = #([1:3;6:8],[2:4;7:9]);","errcatch");
disp(t);

a = []
disp('creation of a as an interval matrix');
ierr($+1) = execstr("a(5,5) = #(5,6);","errcatch");
disp(a);

disp('test for arithmetic');
disp('==================');
disp('addition, substraction :');

disp('|2,3|+|3,4|');
ierr($+1) = execstr("t = #(2,3)+#(3,4);","errcatch");
disp(t);

disp('Matrix with infinite values');
A = #([1,0;5,7],-%inf)
disp('A+2');
ierr($+1) = execstr("t = A + 2","errcatch");
disp(t);
disp('A+0');
ierr($+1) = execstr("t = A + 0","errcatch");
disp(t);
disp('A-2');
ierr($+1) = execstr("t = A - 2","errcatch");
disp(t);

B = #([-4,%inf;0,3],[-%inf,2;6,%inf])
disp('A+B');
ierr($+1) = execstr("t = A + B","errcatch");
disp(t);
disp('A-B');
ierr($+1) = execstr("t = A - B","errcatch");
disp(t);
disp('A+2*B');
ierr($+1) = execstr("t = A + 2* B","errcatch");
disp(t);

disp('==================');
disp('multiplication :');

C = #([-2,-3,2;4,5,0],[-1,2,0;5,6,0])
D = #([0,2;-5,-8;1,10],[0,3;-2,8;2,10.5])
E = D.inf
disp('C * 2');
ierr($+1) = execstr("t = C * 2","errcatch");
disp(t);
disp('C * E');
ierr($+1) = execstr("t = C * E ","errcatch");
disp(t);
disp('C *D');
ierr($+1) = execstr("t = C * D","errcatch");
disp(t);

G = #([-1,-1;0,1],[1,1;0,1])
disp('G * G');
ierr($+1) = execstr("t = G * G","errcatch");
disp(t);
disp('G * G.inf');
ierr($+1) = execstr("t = G * G.inf","errcatch");
disp(t);
disp('G.inf * G');
ierr($+1) = execstr("t = G.inf * G","errcatch");
disp(t);

disp('Matrix with infinite values');

C
disp("");
F = [%inf,12;-10,0;-%inf,%inf]
disp( "C * F");
ierr($+1) = execstr("t = C * F","errcatch");
disp(t);
disp(" [ |-inf,-2| , |-1,inf| ] * |-2,0| = ");
ierr($+1) = execstr("t = #([-%inf,-1],[-2,%inf])*#(-2,0)","errcatch");
disp(t);
disp("");
A
disp("");
B
disp('A * B');
ierr($+1) = execstr("t = A * B","errcatch");
disp(t);
disp('A .* B');
ierr($+1) = execstr("t = A .* B","errcatch");
disp(t);
disp('A .*. B');
ierr($+1) = execstr("t = A .*. B","errcatch");
disp(t);

disp('==================');
disp('division :');

C
disp("");
disp("1 / C")
ierr($+1) = execstr("t = 1 / C ","errcatch");
disp(t);

F = #(F)
disp("");
disp("1 / F")
ierr($+1) = execstr("t = 1 / F ","errcatch");
disp(t);

F = F'
disp("");
disp("C ./ F")
ierr($+1) = execstr("t = C ./ F ","errcatch");
disp(t);

disp("|-4,4| / |-3,9|")

ierr($+1) = execstr("t = #(-4,4) / #(-3,9) ","errcatch");
disp(t);
ierr($+1) = execstr("t = extdiv(#(-4,4) , #(-3,9)) ","errcatch");
disp(t);

disp("|-%inf,4| / |-3,-2|")

ierr($+1) = execstr("t = #(-%inf,4) / #(-3,-2) ","errcatch");
disp(t);

disp("|-%inf,-4| / |-3,2|")

ierr($+1) = execstr("t = #(-%inf,-4) / #(-3,2) ","errcatch");
disp(t);
ierr($+1) = execstr("t = extdiv(#(-%inf,-4),#(-3,2)) ","errcatch");
disp(t);

disp('==================');
disp('roots and power :');

G
disp("");
GG = (G^2 == G*G);
gg = find(~GG);
if gg == []
ierr($+1) = 0;
disp("G^2");
G^2
else
ierr($+1) = 1;
end

disp("");
disp("G^58");
G^58
H = G^58;
if ( H(1,2).inf == -58 & H(1,2).sup == 58 )
ierr($+1) = 0;
else
ierr($+1) = 1;
end

disp("");
A 
disp(" A .^ [1,2;3,4]")
ierr($+1) = execstr("t = A .^ [1,2;3,4]","errcatch");
disp(t);

disp(" A .^ [1,2;1./3,-3./4]")
ierr($+1) = execstr("t = A .^ [1,2;1./3,-3./4]","errcatch");
disp(t);

disp("");
A(2,2) = #(2,3);
A
disp("");
disp(" A .^ [1,2;1./3,-3./4]")
ierr($+1) = execstr("t = A .^ [1,2;1./3,-3./4]","errcatch");
disp(t);

disp('================');
disp('arithmetic functions:');

A = [1,2;-3,-4];

B = [2,3;2,3];
disp("");

C = interval(A,B)

disp('absolute value');
ierr($+1) = execstr("t = abs(C)","errcatch"); //absolute value
disp(t);

disp('exp');
ierr($+1) = execstr("t =  exp(C)","errcatch"); //exponential
disp(t);

disp('log');
ierr($+1) = execstr("t = log(C)","errcatch"); //logarithm
disp(t);

disp('log10');
ierr($+1) = execstr("t = log10(C)","errcatch"); // 10th based logarithm
disp(t);

A = [%pi/6;-%pi/3];

B = [%pi/3;-%pi/6];
disp("");

C = interval(A,B)

disp('sin');
ierr($+1) = execstr("t = sin(C)","errcatch");
disp(t);

disp('cos');
ierr($+1) = execstr("t = cos(C)","errcatch");
disp(t);

disp('tan');
ierr($+1) = execstr("t = tan(C)","errcatch");
disp(t);

disp('cotg');
ierr($+1) = execstr("t = cotg(C)","errcatch");
disp(t);

A = [0.2;0.5];

B = [1;0.9];
disp("");

C = interval(A,B)

disp('asin');
ierr($+1) = execstr("t = asin(C)","errcatch");
disp(t);

disp('acos');
ierr($+1) = execstr("t = acos(C)","errcatch");
disp(t);

disp('atan');
ierr($+1) = execstr("t = atan(C)","errcatch");
disp(t);

A = [2;3];

B = [4;6];
disp("");

C = interval(A,B)

disp('sinh');
ierr($+1) = execstr("t = sinh(C)","errcatch");
disp(t);

disp('cosh');
ierr($+1) = execstr("t = cosh(C)","errcatch");
disp(t);

disp('tanh');
ierr($+1) = execstr("t = tanh(C)","errcatch");
disp(t);

disp('coth');
ierr($+1) = execstr("t = coth(C)","errcatch");
disp(t);

disp('asinh');
ierr($+1) = execstr("t = asinh(C)","errcatch");
disp(t);

disp('acosh');
ierr($+1) = execstr("t = acosh(C)","errcatch");
disp(t);

disp('atanh');
ierr($+1) = execstr("t = atanh(1/C)","errcatch");
disp(t);

disp('==================');
disp('Linear system resolution of A . X = b :');

Ai = [2,0;1,2];
As = [3,1;2,3];
A = interval(Ai,As)
disp('");
bi = [0;60];
bs = [120;240];
b = interval(bi,bs)
disp('");

disp('Default execution ');  
ierr($+1) = execstr("t = A \ b ","errcatch");
disp(t);

disp('Gauss elimination '); 
ierr($+1) = execstr("x = I4Slinearsolve(A,b,""GE"") ","errcatch");
disp(x);

disp(' Try to see if x is the sharpest solution  '); 
ierr($+1) = execstr(" t = I4Slinearsolve(A,b,x,""K"") ","errcatch");
disp(t);

disp ('Other contracted method');
ierr($+1) = execstr(" t = I4Slinearsolve(A,b,x,""GS"")","errcatch");
disp(t);

disp('Let''s try with preconditioning methods'); 
ierr($+1) = execstr("x = I4Slinearsolve(A,b,""PGE"")","errcatch");
disp(x);
disp("");
 
A = #([4,-1,1.5;-0.5,-7,1;-1.5,-0.7,2],[5,1,2.5;.5,-5,2;-.5,-.5,3])
disp('");
b = #([3;0;3],[4;2;4])
disp('");
  
disp('Gauss elimination ');  
ierr($+1) = execstr("t = I4Slinearsolve(A,b,""GE"")","errcatch");
disp(t);

disp('Preconditioned Hansen Bliek');  
ierr($+1) = execstr("t = I4Slinearsolve(A,b,""PHB"")","errcatch");
disp(t);

disp('Default execution ');  
ierr($+1) = execstr("x = A\b","errcatch");
disp(x);

disp('Gauss Seidel ');  
ierr($+1) =  execstr("t = I4Slinearsolve(A,b,x,""GS"")","errcatch");
disp(t);

disp('Gauss Seidel with differents iterative conditions ');  
ierr($+1) =  execstr("t = I4Slinearsolve(A,b,x,""GS"",10,0.89)","errcatch");
disp(t);



q = timer();

disp('==================');
T = (ierr == 0);
ii  = find(~T);

if ii == []
disp (" No problems for Int4Sci !!")
if (q > 1)
disp( " But it seems to be slow !!")
end
else
disp ('There is a problem somewhere !');
printf('The %f th test has failed', ii(1));
disp ('Contact technical support !');
disp('==================');
abort
end 
disp('==================');