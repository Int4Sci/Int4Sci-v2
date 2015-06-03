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

function intervaloperator()
	
	disp(' Interval Arithmetic library is loaded ');
	
endfunction

function intervalLSS()
	
	disp(' Interval Linear System Solver library is loaded ');
	
endfunction

function I4Stest()

	exec(CURRENT_PATH+"/tests/I4Stest.sce");

endfunction

function x = setround(x)

  x = sciroundmode(x);
  
endfunction

//Definition of the type interval Matrix  with inferior and superior bound matrices

function x=interval(varargin)

	
	if length(varargin) == 2
		a=varargin(1);
		b=varargin(2);
	elseif length(varargin)==1
		a=varargin(1);
		b=varargin(1);
	else 
		disp('Too many arguments for creating an interval');
		error('Too many arguments for creating an interval',10001);
		
	end
	
	if typeof(a) == "interM"
	
		x = a;

	else
	
	if size(a) == size(b)

		d  = 0;
		for i = 1:length(a)

			if a(i) >  b(i) 
			
			if  %I4S_arith == "classical"

				c = a(i);
				a(i) = b(i);
				b(i) = c;

			else 
				d=d+1;
		
			end
			end
		end

		if d == 0
			
			x=mlist(["interM","dims","inf","sup","arith","mode"],size(a),a,b,"classical","double");

		else 
	
			x=mlist(["interM","dims","inf","sup","arith","mode"],size(a),a,b,"generalized","double");
		end
				
	else
		if prod(size(a))==1
			x = interval(a*ones(b),b);
		elseif prod(size(b))==1
			x = interval(a,b*ones(a));
		else
			disp('Invalid size for arguments');
			error('Invalid size for arguments',10002);
		end
	end
	end

endfunction

function y=#(varargin)

	if length(varargin) == 2
		a=varargin(1);
		b=varargin(2);
	elseif length(varargin)==1
		a=varargin(1);
		b=varargin(1);
	else 
		disp('Too many arguments for creating an interval');
		error('Too many arguments for creating an interval',10001);
	end

	y=interval(a,b);

endfunction

// check the size of an interval Matrix

function R=%interM_size(varargin)

	if argn(2) == 2
		M=varargin($-1);
		flag=varargin($);
		
		if flag == 1 | flag == "r"
			R=M.dims(1);
		elseif flag == 2 | flag == "c"
			R=M.dims(2);
		elseif flag == "*"
			R=M.dims(1)*M.dims(2);
		else
			disp('Incorrect function arguments');
			error('Incorrect function arguments', 10004);
		end
	
	elseif argn(2) == 1
		R=varargin(1).dims;

	else
		disp('Incorrect function arguments');
		error('Incorrect function arguments', 10004);	
	end
	
endfunction

//Definition of the type interval Matrix with mid value and radium

function xe=mradius(x,e)
	
	if ~exists('x') | ~exists('e')
		error(39);
	else
	a = size(x);
	b = size(e);
	c = [1,1];

	if (a(1) <> b(1) | a(2) <> b(2)) & prod(b) <> 1
		disp('Wrong radius size');
		error('Wrong radius size',10003);
	else

  		if typeof(x) == 'interM'
    			xe=interval(x.inf-abs(e),x.sup+abs(e));
 	 	else
    			xe=interval(x-abs(e),x+abs(e));
 	 	end
	end
	end

endfunction


//Display print of a interval Matrix

//Automaticly
function %interM_p(x)
	
	if prod(x.dims) > 100000
		disp ('Sorry, the printing is too long !!!');		

	else

	if x.arith == "generalized"
		if %I4S_arith == "generalized"
			disp('Warning : Generalized interval');
		else
			x.arith = "classical";
		end
	end	
	
	st = "|" + string(x.inf) + "," + string(x.sup) + "|";

	st= matrix(st,x.dims);
	disp(st);

	end

endfunction

//Via string function
function st= %interM_string(x)
	
	st = "|" + string(x.inf) + "," + string(x.sup) + "|";

	st = matrix(st,x.dims);


endfunction

// extraction from an interval matrix

function M=%interM_e(varargin)

  	M=varargin($);
  	ind=matrix(1:prod(M.dims),M.dims);
  	ind=ind(varargin(1:$-1));

  	if size(ind,"*") ==0 then
   		 M=[];
 	 else
   		 M=interval(matrix(M.inf(ind),size(ind)),matrix(M.sup(ind),size(ind)))
 	 end

endfunction


// column concatenation

function R=%interM_c_interM(varargin)
 	 R=interval([varargin($-1).inf,varargin($).inf],[varargin($-1).sup,varargin($).sup]); 
endfunction

// row concatenation

function R=%interM_f_interM(varargin)

  	R=interval([varargin($-1).inf;varargin($).inf],[varargin($-1).sup;varargin($).sup]);
 
endfunction


// column concatenation

function R=%interM_c_s(varargin)

  	R=interval([varargin($-1).inf,varargin($)],[varargin($-1).sup,varargin($)]); 

endfunction

// row concatenation

function R=%interM_f_s(varargin)

 	 R=interval([varargin($-1).inf;varargin($)],[varargin($-1).sup;varargin($)]); 

endfunction


// column concatenation

function R=%s_c_interM(varargin)

  	R=interval([varargin($-1),varargin($).inf],[varargin($-1),varargin($).sup]); 

endfunction

// row concatenation

function R=%s_f_interM(varargin)

  	R=interval([varargin($-1);varargin($).inf],[varargin($-1);varargin($).sup]); 

endfunction

//Insertion in an interval matrix

//Interval submatrix insertion

function R=%interM_i_interM(varargin)  
  
	M = varargin($);                           
	M.inf(varargin(1:$-2)) = varargin($-1).inf;
	M.sup(varargin(1:$-2)) = varargin($-1).sup;
	R = #(M.inf,M.sup);                       

endfunction


//Real submatrix insertion

function R=%s_i_interM(varargin)
	
	M=varargin($);
  	x=varargin($-1);

	if prod(size(x)) == 0 & prod(size(M(varargin(1:$-2)))) == prod(size(M))

        	R=interval(%nan);	
    	else
		M.inf(varargin(1:$-2)) = x;
		M.sup(varargin(1:$-2)) = x;
		R = #(M.inf,M.sup);                       
	
  	 end 

endfunction

//Insertion of interval submatrix in an real matrix
function M=%interM_i_s(varargin)

	if varargin($) == []

        	M=interval(0);
    	else

        	M=interval(varargin($));

   	 end 

  	R=varargin($-1);

  	if length(varargin) == 4

   		 M(varargin($-3),varargin($-2)) = R;  

  	elseif length(varargin) == 3

    		M(varargin($-2)) = R;

 	 end

endfunction

//transpose

function R=%interM_t(M)

	R=interval(M.inf',M.sup');

endfunction

//length


function R=lengthM(varargin)
	if typeof(varargin(1)) == 'interM' then

		if argn(2) == 1 then
		
			M = varargin(1);
			R = prod(M.dims);

		else

			disp('The inputs of the function are incorrect ');
			error('The inputs of the function are incorrect ', 10004);

		end
	
	else

		R = length(varargin(1:$));
	end
endfunction

//function s=oldlength(e)
//	s=length(e);
//endfunction

//bisection

function R=ncut_bisect(d,val,nu)

       t=d(:);
       R=list();

   	if val >= 1 & (nu <= prod(d.dims) & nu >0) & length(val) == 1

      		b=t(nu).inf:width(t(nu))/val:t(nu).sup;
		if int(val)==val	
		   b($)=t(nu).sup;
		elseif  b($) ~=  t(nu).sup  
			b($+1) = t(nu).sup; 
      		end 
       		if b($) ~=  t(nu).sup
              		b($+1) = t(nu).sup;
      		 end 

         	 for i=1:length(b)-1               
			t=d;
          		t.inf(nu)=b(i);
           		t.sup(nu)=b(i+1);
           		R(i)=interval(matrix(t.inf,d.dims),matrix(t.sup,d.dims));
         	end

   	elseif (val >0 & val < 1 ) & (nu <= prod(d.dims) & nu >0) & length(val) == 1

       		t=d;
       		t.sup(nu)=t.inf(nu)+width(t(nu))*val;
       		R(1)=interval(matrix(t.inf,d.dims),matrix(t.sup,d.dims));
       		t=d;
       		t.inf(nu)=t.inf(nu)+width(t(nu))*val;
       		R(2)=interval(matrix(t.inf,d.dims),matrix(t.sup,d.dims));
       
	else
       	R=list();
	disp('Incorrect arguments for bisection');
      	error('Incorrect arguments for bisection', 10006);

   	end

endfunction 


function R=ncut_bisect_list(L,val,nu)

	R=ncut_bisect(L(1),val,nu);
   	for i=2:lengthM(L)
      	 T=ncut_bisect(L(i),val,nu);
      	 for j=1:lengthM(T)
           R($+1)=T(j);
       	end
   	end

endfunction

function R=%interM_b_s(varargin)

   	t=varargin(1);
  	c=varargin(2);


  	 if length(varargin) ==2

      	 	[val,nu]=max(width(t(:))); // cut on largest dimension
		
		tt = t;
		ii = 0;
		 while isinf(val)
			ii = ii + 1;
			if ii > lengthM(tt)-1;	
				disp ('function/operation does not handle infinite value');
				error ('function/operation does not handle infinite value', 20008);
			end
			tt(nu) = #(0);
			[val,nu] = max(width(tt(:)));	
		end
		
      	 	R=ncut_bisect(t,c,nu);

  	elseif length(varargin) ==3

       		nnu=varargin(3);
	
		if isinf(width(t(nnu)))
			disp ('function/operation does not handle infinite value');
			error ('function/operation does not handle infinite value', 20008);
		end
		
		 if size(nnu,'r') ==1 | size(nnu,'c') ==1  
                   
	 	R=ncut_bisect(t,c,nnu(1));
           	for i=2:lengthM(nnu)
              	 R=ncut_bisect_list(R,c,nnu(i));    
          	end
                     
       		else
          	 R=list();
		disp('Incorrect arguments for bisection');
           	error('Incorrect arguments for bisection', 10006);

       		end

   	else
		disp('Incorrect arguments for bisection');
       		error('Incorrect arguments for bisection', 10006);
   	end

endfunction 

////////////////////////////////////////////////////
// tools

// order interv double

function B=%interM_1_s(varargin)

  B = varargin(1).sup < varargin($);

endfunction

function B=%interM_3_s(varargin)

  B = varargin(1).sup <= varargin($);

endfunction

function B=%interM_2_s(varargin)

  B = varargin(1).inf > varargin($);

endfunction

function B=%interM_4_s(varargin)

  B = varargin(1).inf >= varargin($);

endfunction

// order double interval

function B=%s_1_interM(varargin)

  B = varargin(1) < varargin($).sup;

endfunction

function B=%s_3_interM(varargin)

  B = varargin(1) <= varargin($).sup;

endfunction

function B=%s_2_interM(varargin)

  B = varargin(1) > varargin($).inf;

endfunction

function B=%s_4_interM(varargin)

  B = varargin(1) >= varargin($).inf;

endfunction

//order interval - interval

function B=%interM_1_interM(varargin)

  B = varargin(1).sup < varargin($).inf;

endfunction

function B=%interM_3_interM(varargin)

  B = varargin(1).sup <= varargin($).inf;

endfunction

function B=%interM_2_interM(varargin)

  B = varargin(1).inf > varargin($).sup;

endfunction

function B=%interM_4_interM(varargin)

  B = varargin(1).inf >= varargin($).sup;

endfunction

// == & <>

function B=%interM_o_interM(a,b)

	 B = (a.inf == b.inf & a.sup == b.sup);
	
endfunction

function B=%interM_n_interM(a,b)

	 B = (a.inf ~= b.inf | a.sup ~= b.sup);
	
endfunction	

function B=%interM_o_s(a,b)

	 B = (a.inf == b & a.sup == b);
	
endfunction

function B=%interM_n_s(a,b)

	 B = (a.inf ~= b | a.sup ~= b);
	
endfunction	

function B=%s_o_interM(b,a)

	 B = (a.inf == b & a.sup == b);
	
endfunction

function B=%s_n_interM(b,a)

	 B = (a.inf ~= b | a.sup ~= b);
	
endfunction	

//abs - magnitude

function R=mag(varargin)

  R=max(abs(varargin(1).inf),abs(varargin(1).sup))

endfunction

// mignitude

function R=mig(varargin)

  for i=1:prod(size(varargin(1)))

    if varargin(1).inf(i) > 0 then

      R(i)=varargin(1).inf(i)

    elseif varargin(1).sup(i) < 0 then

      R(i)=abs(-varargin(1).sup(i))

    else

      R(i)=0;

    end

  end 
 
  R=matrix(R,size(varargin(1)));

endfunction

//norm

function R=%interM_norm(varargin)

  if argn(2) == 2 then
	if length(varargin(1).dims)==2 & min(varargin(1).dims)==1 then
		
		if varargin(2) == 'inf' |  varargin(2) == %inf
	
			R = max(mag(varargin(1)));

		elseif typeof(varargin(2)) == "constant" | length(varargin(2)) == 1 then

			p = varargin(2);	
			R = abs(sum(varargin(1)^p))^(1./p);

		else
			R=sqrt(sum(varargin(1)^2));

 		 end

	else
		if varargin(2) == 'inf' |  varargin(2) == %inf then
			
			for i=1:size(varargin(1),'r')

				D(i)=sum(mag(varargin(1)(i,:)));
			end
		
			R=max(D);

		elseif varargin(2) == 1 then
	
			for i=1:size(varargin(1),'c')

				D(i)=sum(mag(varargin(1)(:,i)));
			end
		
			R=max(D);

		else

			for i=1:size(varargin(1),'r')

				D(i)=sum(mag(varargin(1)(i,:)));
			end
		
			R=max(D);	

 		 end

	end
		
   else

	if length(varargin(1).dims)==2 & min(varargin(1).dims)==1 then

		norm(varargin(1),2);
	else

		norm(varargin(1),'inf');
	end
  end

endfunction

function R=%interM_sum(varargin)

  if argn(2) == 2 then

	 if varargin(2) == 'r' | varargin(2) == 1 then
	
		for  i=1:size(varargin(1),'c')
		
			R(i) = sum(varargin(1)(:,i));
		end

		R = matrix(R,[1,size(varargin(1),'c')]);
	
   	elseif varargin(2) == 'c' | varargin(2) == 2 then

		for  i=1:size(varargin(1),'r')
		
			R(i) = sum(varargin(1)(i,:));
		end
 
   	elseif varargin(2) == 'm'

		if varargin(1).dims(1) == 1

			R = sum(varargin(1),'c');
		else
			if  varargin(1).dims(2) != 1

				R = sum(varargin(1),'r');
			else

				R = varargin(1);
			end
		end
  
   	 else

    		R=interval(sum(varargin(1).inf),sum(varargin(1).sup));

    	end

  else

	R=interval(sum(varargin(1).inf),sum(varargin(1).sup));

  end
 
endfunction

//diagonal

function R=%interM_diag(varargin)

  if argn(2) == 2 then

	R = interval(diag(varargin(1).inf,varargin(2)),diag(varargin(1).sup,varargin(2)));

  else

	R = interval(diag(varargin(1).inf),diag(varargin(1).sup));	
  end

endfunction

// middle of an inteval

function m=mid(x)

   m=(x.sup+x.inf)./2;

endfunction

// radius of an interval

function m=rad(x)

   m=0.5*abs(x.sup-x.inf)

endfunction

//width of an interval

function m=width(x)

   m=abs(x.sup-x.inf)

endfunction

function m=wid(x)

   m=width(x)

endfunction

//supbound & infbound of an interval

function m=sup(x)

	m=x.sup

endfunction

function m=inf(x)

	m=x.inf

endfunction

function c=intersection(a,b)


  if typeof(a) == 'interM'
	 if a.arith == "generalized" 

		NIY();
	end
  end	

  if typeof(b) == 'interM' 
	if b.arith == "generalized" 

		NIY();
	end
  end	

  if size(a) == size(b)

    if typeof(a)<> 'interM'

      a=interval(a);

    end

    if typeof(b)<> 'interM'

      b=interval(b);

    end

    c=a;

    for i=1:prod(a.dims)
      
	binf=max(a.inf(i),b.inf(i));
	bsup=min(a.sup(i),b.sup(i));
	
	if binf > bsup
	  
  		c(i)=interval(%nan);

	else

	  c(i)=interval(binf,bsup);

	end
	
    end

  else 

    disp('dimension problem !! try  with two similar dimension intervals');
    error('impossible operation',20005);

  end

endfunction

function c=intersec(a,b)
    
    c=intersection(a,b)

endfunction

function c = hull(varargin)

	if lengthM(varargin) == 1
		
		a = varargin(1);

		if typeof(a) == "list"
		
		b=[];
		for i = 1:lengthM(a)

			ierr = execstr('b = [b;a(i)];','errcatch');
			if ierr <> 0
				disp('Input dimensions problem !! ');
				error('impossible operation',20005);	
			end		
		end
		a = b;
		end

		if prod(size(a)) == 1

			c = #(a);

		else
			
			b = #(a(:));
			c = #(min(inf(b)),max(sup(b)));

		end

	else 
		if length(varargin) == 2

		a = varargin(1);
		b = varargin(2);
		
		if typeof(a) == "list"
		
		f=[];
		for i = 1:lengthM(a)

			ierr = execstr('f = [f;a(i)];','errcatch');
			if ierr <> 0
				disp('Input dimensions problem !! ');
				error('impossible operation',20005);	
			end		
		end
		a = f;
		end

		if typeof(b) == "list"
		
		g=[];
		for i = 1:lengthM(b)

			ierr = execstr('g = [g;b(i)];','errcatch');
			if ierr <> 0
				disp('Input dimensions problem !! ');
				error('impossible operation',20005);	
			end		
		end
		b = g;
		end

		e = [inf(#(a(:)));inf(#(b(:)))];
		d = [sup(#(a(:)));sup(#(b(:)))];

		c = #(min(e),max(d));

		else

		error(39);
		
		end

	end
		

endfunction

function c = dothull(a,b)

	if typeof(a) == "list"
		
		f=[];
		for i = 1:lengthM(a)

			ierr = execstr('f = [f;a(i)];','errcatch');
			if ierr <> 0
				disp('Input dimensions problem !! ');
				error('impossible operation',20005);	
			end		
		end
		a = f;
	end

	if typeof(b) == "list"
		
		g=[];
		for i = 1:lengthM(b)

			ierr = execstr('g = [g;b(i)];','errcatch');
			if ierr <> 0
				disp('Input dimensions problem !! ');
				error('impossible operation',20005);	
			end		
		end
		b = g;
	end

	if lengthM(a) == lengthM(b) & size(a,'r') == size(b,'r')

		for i = 1:lengthM(a)
			c(i) = hull(a(i),b(i));
		end

		c = matrix(c,size(a));
	
	elseif lengthM(a) == 1
		
		for i = 1:lengthM(b)
			c(i) = hull(a,b(i));
		end

		c = matrix(c,size(b));	

	elseif lengthM(b) == 1
			
		for i = 1:lengthM(a)
			c(i) = hull(a(i),b);
		end

		c = matrix(c,size(a));	

	else
		disp('Input dimensions problem !! ');
		error('impossible operation',20005);	
	end
endfunction

function r =isemptyset(a)

if typeof(a) == 'interM'	
	rinf = isnan(a.inf) ;
	rsup = isnan(a.sup) ;

	if isnan(a.inf) == isnan(a.sup)

		r = rinf;

	else 

		for i = 1:prod(a.dims)
			if rinf(i) ~= rsup(i)
			ri(i) = %F;
			else	
			ri(i) = rinf(i);
			end	
		end

		r = matrix(ri,a.dims);

	end
else

	r = isempty(a);

end

endfunction


function y=%interM_matrix(varargin)

  y=interval(matrix(varargin(1).inf,varargin(2:$)),matrix(varargin(1).sup,varargin(2:$)));  

endfunction

// return a boolean matrix to check if x is inside y

function res=in(x,y)

  if (size(x,'r') <> size(y,'r') | size(x,'c') <> size(y,'c')) & size(x,"*") <> 1

    disp('Incorrect matrix size in in/strictin function');	
    error("Incorrect matrix size in in/strictin function",10010);

  end

  if typeof(y) <> 'interM'

    disp('Not defined inside the function IN');
    error("Not defined inside the function IN",20011);

  end

  if typeof(x) == 'interM'

    res=matrix((x.sup(:) <= y.sup(:) & x.inf(:)>=y.inf(:)),x.dims);

  elseif typeof(x) == 'constant'

    res=matrix((x(:) <= y.sup(:) & x(:)>=y.inf(:)),y.dims);

  else
    
    disp('Not defined inside the function IN');
    error("Not defined inside the function IN",20011);
  end

endfunction

// return a boolean matrix to chech if x is strictly inside y

function res=strictin(x,y)

  if (size(x,'r') <> size(y,'r') | size(x,'c') <> size(y,'c')) & size(x,"*") <> 1

    disp('Incorrect matrix size in in/strictin function');	
    error("Incorrect matrix size in in/strictin function",10010);

  end

  if typeof(y) <> 'interM'

    disp('Not defined inside the function IN');
    error("Not defined inside the function IN",20011);

  end

  if typeof(x) == 'interM'

    res=matrix((x.sup(:) < y.sup(:) & x.inf(:) > y.inf(:)),x.dims);

  elseif typeof(x) == 'constant'

    res=matrix((x(:) < y.sup(:) & x(:) > y.inf(:)),y.dims);

  end

endfunction

//////////////////////////////////
//Basic Errors handling

//Not Implented Yet Error message

function NIY()

	disp('Operation not yet implemented');
	error('Operation not yet implemented',100000);

endfunction



//Internal Error Message

function IEM()

	disp('Internal problem ! Contact technical support !');
	error('Internal Problem',90000);

endfunction

// Flag reading

function FlagReading(flag)

	if (flag == -1)
	
		disp('Input dimensions problem !! ');
		error('impossible operation',20005);

	end

	if flag == -2

		IEM();

	end

	if flag == -3

		NIY();

	end
	
	if flag == -4
	
		disp('Allocation Problem !');
		disp('Free memory needed !! ');
		error('Allocation problem',90014);
	end

	if (flag == -5)

 		disp('Problem when solving a linear system');
		disp('Maybe this system has no solution or a singularity !');
		disp('Try an other method and see help files');
		error('Problem when solving a linear system',30015);

	end

	if (flag == -6)

		disp('Incorrect dimensions for a linear system');
		error('Incorrect dimensions for a linear system',30016);

	end

	if (flag == -7)

		disp('Hansen Bliek method doesn''t work !');
		disp('Try an other method and see help files');
		error('Hansen Bliek method doesn''t work !',10020);

	end

endfunction

function [aa,bb] = ArithReading(a,b)

  c = 0;
  d = 0;
 
	if typeof(a) == "constant"

		a = interval(a);
    		c = 1;
    
	end

	if typeof(b) == "constant"

		b = interval(b);
		d = 1;

	end

	if typeof(a) ~= "interM" | typeof(b) ~= "interM"
	
		disp('Input types problem !! ');
		error('impossible operation',20005);	
	end

	if size(a.inf,'r') ~= size(a.sup,'r') | size(a.inf,'c') ~= size(a.sup,'c') 

		if prod(size(a.inf)) == 1 | prod(size(a.sup)) == 1 

			a = interval(a.inf,a.sup);
		
		else
			disp('Input dimensions problem between bounds !! ');
			error('impossible operation',20005);
		
		end
	end
			
	if size(b.inf,'r') ~= size(b.sup,'r') | size(b.inf,'c') ~= size(b.sup,'c') 

		if prod(size(b.inf)) == 1 | prod(size(b.sup)) == 1 

			b = interval(b.inf,b.sup);
		
		else
			disp('Input dimensions problem between bounds !! ');
			error('impossible operation',20005);
		
		end
	end

  x = a.inf > a.sup;
  ii = find(x);
  
  y = b.inf > b.sup;
  jj = find(y);
  
	if a.arith == "generalized" | b.arith == "generalized" | ~(ii == []) | ~(jj == [])
	
		if %I4S_arith == "generalized"

			NIY();
		else
	     
			disp("Warning : redefining your interval for classical arithmetic"); 

			if a.arith == "generalized" | ~(ii == [])
			
				a = interval(a.inf,a.sup);

			end

			if b.arith == "generalized" | ~(jj == [])
			
				b = interval(b.inf,b.sup);

			end

		end

  	end
  	
  	if c == 1
  	
  	   a = a.inf;
  	   
  	end
  	
  	if d == 1
  	
  	 b = b.inf;
  	 
  	end

  aa = a;
  bb = b;
  
endfunction


//////////////////////////////////////
//arithmetique:
//-------------

function c=%interM_s(a)

 c = 0 - a;

endfunction

//Addition of two intervals

function c=%interM_a_interM(a,b)
 		
	[a,b] = ArithReading(a,b);

	[resinf,ressup,flag] = scioperator(a.inf,a.sup,'a',b.inf,b.sup);

	FlagReading(flag);

	if (prod(size(resinf))==0)

		c=[];

	else

		c=interval(resinf,ressup);

	end	

endfunction

//Addition of an interval and a constant

function c=%interM_a_s(a,b)

	[a,b] = ArithReading(a,b);

	[resinf,ressup,flag]= scioperator(a.inf,a.sup,'a',b,[]);

	FlagReading(flag);

	if (prod(size(resinf))==0)

		c=[];

	else

		c=interval(resinf,ressup);

	end

endfunction

function c=%s_a_interM(b,a)

	[a,b] = ArithReading(a,b);

	[resinf,ressup,flag]=scioperator(a.inf,a.sup,'a',b,[]);

	FlagReading(flag);

	if (prod(size(resinf))==0)

		c=[];

	else

		c=interval(resinf,ressup);

	end

endfunction

//Substraction of two intervals

function c=%interM_s_interM(a,b)
 	
	[a,b] = ArithReading(a,b);

	[resinf,ressup,flag] = scioperator(a.inf,a.sup,'s',b.inf,b.sup);

	FlagReading(flag);

	if (prod(size(resinf))==0)

		c=[];

	else

		c=interval(resinf,ressup);

	end	

endfunction

//Substraction of an interval and a constant

function c=%interM_s_s(a,b)

	[a,b] = ArithReading(a,b);

	[resinf,ressup,flag]= scioperator(a.inf,a.sup,'s',b,[]);

	FlagReading(flag);

	if (prod(size(resinf))==0)

		c=[];

	else

		c=interval(resinf,ressup);

	end

endfunction

function c=%s_s_interM(b,a)

	[a,b] = ArithReading(a,b);

	[resinf,ressup,flag]=scioperator(a.inf,a.sup,'s',b,[]);

	FlagReading(flag);

	if (prod(size(resinf))==0)

		c=[];

	else

		c=interval(-ressup,-resinf);

	end

endfunction

// Interval Matrix multiplication

// Multiplication of two interval matrices

function c=%interM_m_interM(a,b)
 	
	[a,b] = ArithReading(a,b);

	[resinf,ressup,flag] = scioperator(a.inf,a.sup,'m',b.inf,b.sup);

	FlagReading(flag);
	
	if (prod(size(resinf))==0)

		c=[];

	else

		c=interval(resinf,ressup);

	end	

endfunction

// Multiplication of an interval matrix and a real matrix

function c=%interM_m_s(a,b)

	[a,b] = ArithReading(a,b);

	[resinf,ressup,flag] = scioperator(a.inf,a.sup,'m',b,[]);

	FlagReading(flag);

	if (prod(size(resinf))==0)

		c=[];

	else

		c=interval(resinf,ressup);

	end

endfunction

function c=%s_m_interM(a,b)

	[a,b] = ArithReading(a,b);

	[resinf,ressup,flag] = scioperator(a,[],'m',b.inf,b.sup);

	FlagReading(flag);

	if (prod(size(resinf))==0)

		c=[];

	else

		c=interval(resinf,ressup);

	end

endfunction	

// Element by element multiplication of two interval matrices

function c=%interM_x_interM(a,b)
 	
	[a,b] = ArithReading(a,b);

	[resinf,ressup,flag] = scioperator(a.inf,a.sup,'x',b.inf,b.sup);

	FlagReading(flag);

	if (prod(size(resinf))==0)

		c=[];

	else

		c=interval(resinf,ressup);

	end	

endfunction

// Element by element multiplication of matrix interval and real matrix

function c=%interM_x_s(a,b)
 	
	[a,b] = ArithReading(a,b);

	[resinf,ressup,flag] = scioperator(a.inf,a.sup,'x',b,[]);

	FlagReading(flag);

	if (prod(size(resinf))==0)

		c=[];

	else

		c=interval(resinf,ressup);

	end	

endfunction

function c=%s_x_interM(a,b)
 	
	[a,b] = ArithReading(a,b);

	[resinf,ressup,flag] = scioperator(b.inf,b.sup,'x',a,[]);

	FlagReading(flag);

	if (prod(size(resinf))==0)

		c=[];

	else

		c=interval(resinf,ressup);

	end	

endfunction


// Kronecker Multiplication of two interval matrices

function c=%interM_k_interM(a,b)
 	
	[a,b] = ArithReading(a,b);

	[resinf,ressup,flag] = scioperator(a.inf,a.sup,'k',b.inf,b.sup);

	if (prod(size(resinf))==0)

		c=[];

	else

		c=interval(resinf,ressup);

	end	

endfunction

// Kronecker Multiplication of an interval matrix and a real matrix

function c=%interM_k_s(a,b)

	[a,b] = ArithReading(a,b);

	[resinf,ressup,flag] = scioperator(a.inf,a.sup,'k',b,[]);

	if (prod(size(resinf))==0)

		c=[];

	else

		c=interval(resinf,ressup);

	end

endfunction

function c=%s_k_interM(a,b)

	[a,b] = ArithReading(a,b);

	[resinf,ressup,flag] = scioperator(a,[],'k',b.inf,b.sup);

	if (prod(size(resinf))==0)

		c=[];

	else

		c=interval(resinf,ressup);

	end

endfunction

// Interval Matrix division

// Division of two interval matrices

function c=%interM_r_interM(a,b)

	[a,b] = ArithReading(a,b);

	if lengthM(a) > lengthM(b) & lengthM(b) == 1 & b == interval(0)

		c = ones(a.inf);
		c(:) = #(%nan,%nan);
		return;
	
	end

	ii = [];
	d = 0;

	for i=1:lengthM(b)
	
		if b(i) == interval(0)
		
			ii = [i,ii];
			b(i) = interval(1);
			d=d+1;

		end
	end
		
	[resinf,ressup,flag] = scioperator(a.inf,a.sup,'r',b.inf,b.sup);

	FlagReading(flag);
	
	if (prod(size(resinf))==0)

		c=[];

	else

		c=interval(resinf,ressup);

		if d ~= 0
			c(ii) = #(%nan,%nan);
		end

	end		

endfunction

//Division of interval by real

function c=%interM_r_s(a,b)
 	
	[a,b] = ArithReading(a,b);

	if lengthM(a) > lengthM(b) & lengthM(b) == 1 & b == 0

		c = ones(a.inf);
		c(:) = #(%nan,%nan);
		return;
	
	end

	ii = [];
	d = 0;

	for i=1:lengthM(b)
	
		if b(i) == 0
		
			ii = [i,ii];
			b(i) = 1;
			d=d+1;
		end
	end

	if size(b,1) == size(b,2) & lengthM(a) == lengthM(b) & lengthM(b) ~= 1

		c = a * b^(-1);

	else

		[resinf,ressup,flag] = scioperator(a.inf,a.sup,'r',b,[]);

		FlagReading(flag);
	
		if (prod(size(resinf))==0)

			c=[];

		else

			c = interval(resinf,ressup);

			if d ~= 0
				c(ii) = #(%nan,%nan);
			end

		end
	end	

endfunction

function c=%s_r_interM(a,b)

	[a,b] = ArithReading(a,b);

	if lengthM(a) > lengthM(b) & lengthM(b) == 1 & b == interval(0)

		c = ones(a);
		c(:) = #(%nan,%nan);

		return;
	
	end

	ii = [];
	d = 0;

	for i=1:lengthM(b)
	
		if b(i) == interval(0)
		
			ii = [i,ii];
			b(i) = interval(1);
			d = d+1;

		end
	end

	[resinf,ressup,flag] = scioperator(a,[],'r',b.inf,b.sup);

	FlagReading(flag);
	
	if (prod(size(resinf))==0)

		c=[];

	else

		c=interval(resinf,ressup);

		if d ~= 0
			c(ii) = #(%nan,%nan);

		end

	end	

endfunction

// Left Division of two interval matrices
// or interval linear system resolution

function c=%interM_l_interM(a,b)

	[a,b] = ArithReading(a,b);
	 
  if (size(b,'r') == size(a,'c') & size(b,'c') == 1)
  
    c = I4Slinearsolve(a,b);
  
  else
	if lengthM(a) < lengthM(b) & lengthM(a) == 1 & a == interval(0)

		c = ones(b.inf);
		c(:) = #(%nan,%nan);

		return;
	
	end

	ii = [];
	d = 0;

	for i=1:lengthMM(a)
	
		if a(i) == interval(0)
		
			ii = [i,ii];
			a(i) = interval(1);
			d = d+1;

		end
	end

	[resinf,ressup,flag] = scioperator(b.inf,b.sup,'r',a.inf,a.sup);

	FlagReading(flag);
	
	if (prod(size(resinf))==0)

		c=[];

	else

		c=interval(resinf,ressup);

		if d ~= 0
			c(ii) = #(%nan,%nan);

		end

	end	
	end

endfunction

//Left Division of interval by real
// or interval linear system resolution

function c=%interM_l_s(a,b)
 	
	[a,b] = ArithReading(a,b);

  if (size(b,'r') == size(a,'c') & size(b,'c') == 1)
  
    c = I4Slinearsolve(a,b);
  
  else
	if lengthM(b) > lengthM(a) & lengthM(a) == 1 & a == interval(0)

		c = ones(b);
		c(:) = #(%nan,%nan);

		return;
	
	end

	ii = [];
	d = 0;

	for i=1:lengthM(a)
	
		if a(i) == interval(0)
		
			ii = [i,ii];
			a(i) = interval(1);
			d=d+1;
		end
	end

	if size(a,1) == size(a,2) & lengthM(a) == lengthM(b) & lengthM(b) ~= 1

		c = b * a^(-1);

	else

		[resinf,ressup,flag] = scioperator(b,[],'r',a.inf,a.sup);

		FlagReading(flag);
	
		if (prod(size(resinf))==0)

			c=[];

		else

			c=interval(resinf,ressup);

			if d ~= 0
				c(ii) = #(%nan,%nan);

			end

		end
	end	
	end

endfunction

function c=%s_l_interM(a,b)

	[a,b] = ArithReading(a,b);

  if (size(b,'r') == size(a,'c') & size(b,'c') == 1)
  
    c = I4Slinearsolve(a,b);
  
  else
	if lengthM(b) > lengthM(a) & lengthM(a) == 1 & a == 0

		c = ones(b.inf);
		c(:) = #(%nan,%nan);

		return;
	
	end

	ii = [];
	d = 0;

	for i=1:lengthM(b)
	
		if a(i) == 0
		
			ii = [i,ii];
			a(i) = 1;
			d=d+1;
		end
	end

	[resinf,ressup,flag] = scioperator(b.inf,b.sup,'r',a,[]);

	FlagReading(flag);
	
	if (prod(size(resinf))==0)

		c=[];

	else

		c=interval(resinf,ressup);

		if d ~= 0
			c(ii) = #(%nan,%nan);

		end

	end	
	end

endfunction

// Element by element Interval Matrix division

// Element by element Division of two interval matrices

function c=%interM_d_interM(a,b)
 	
	[a,b] = ArithReading(a,b);

	if lengthM(a) > lengthM(b) & lengthM(b) == 1 & b == interval(0)

		c = ones(a.inf);
		c(:) = #(%nan,%nan);

		return;
	
	end

	ii = [];
	d = 0;

	for i=1:lengthM(b)
	
		if b(i) == interval(0)
		
			ii = [i,ii];
			b(i) = interval(1);
			d=d+1;

		end
	end

	[resinf,ressup,flag] = scioperator(a.inf,a.sup,'d',b.inf,b.sup);

	FlagReading(flag);
	
	if (prod(size(resinf))==0)

		c=[];

	else

		c=interval(resinf,ressup);

		if d ~= 0
			c(ii) = #(%nan,%nan);

		end

	end	

endfunction

// Element by element Division of interval by real

function c=%interM_d_s(a,b)
 	
	[a,b] = ArithReading(a,b);

	if lengthM(a) > lengthM(b) & lengthM(b) == 1 & b == 0

		c = ones(a.inf);
		c(:) = #(%nan,%nan);

		return;
	
	end

	ii = [];
	d = 0;

	for i=1:lengthM(b)
	
		if b(i) == 0
		
			ii = [i,ii];
			b(i) = 1;
			d=d+1;
		end
	end

	[resinf,ressup,flag] = scioperator(a.inf,a.sup,'d',b,[]);

	FlagReading(flag);
	
	if (prod(size(resinf))==0)

		c=[];

	else

		c=interval(resinf,ressup);

		if d ~= 0
			c(ii) = #(%nan,%nan);

		end

	end	

endfunction

function c=%s_d_interM(a,b)

	[a,b] = ArithReading(a,b);

	if lengthM(a) > lengthM(b) & lengthM(b) == 1 & b == interval(0)

		c = ones(a);
		c(:) = #(%nan,%nan);

		return;
	
	end

	ii = [];
	d = 0;

	for i=1:lengthM(b)
	
		if b(i) == interval(0)
		
			ii = [i,ii];
			b(i) = interval(1);
			d = d+1;

		end
	end

	[resinf,ressup,flag] = scioperator(a,[],'d',b.inf,b.sup);

	FlagReading(flag);
	
	if (prod(size(resinf))==0)

		c=[];

	else

		c=interval(resinf,ressup);
	
		if d ~= 0
			c(ii) = #(%nan,%nan);

		end

	end	

endfunction

//  Element by element Left Division of two interval matrices

function c=%interM_q_interM(a,b)

	[a,b] = ArithReading(a,b);

	if lengthM(b) > lengthM(a) & lengthM(a) == 1 & a == interval(0)

		c = ones(b);
		c(:) = #(%nan,%nan);

		return;
	
	end

	ii = [];
	d = 0;

	for i=1:lengthM(a)
	
		if a(i) == interval(0)
		
			ii = [i,ii];
			a(i) = interval(1);
			d=d+1;
		end
	end

	[resinf,ressup,flag] = scioperator(b.inf,b.sup,'d',a.inf,a.sup);

	FlagReading(flag);
	
	if (prod(size(resinf))==0)

		c=[];

	else

		c=interval(resinf,ressup);

		if d ~= 0
			c(ii) = #(%nan,%nan);

		end

	end	

endfunction

//Elementwise Left Division of interval by real

function c=%interM_q_s(a,b)
 	
	[a,b] = ArithReading(a,b);

	if lengthM(b) > lengthM(a) & lengthM(a) == 1 & a == interval(0)

		c = ones(b);
		c(:) = #(%nan,%nan);

		return;
	
	end

	ii = [];
	d = 0;

	for i=1:lengthM(a)
	
		if a(i) == interval(0)
		
			ii = [i,ii];
			a(i) = interval(1);
			d=d+1;
		end
	end

	[resinf,ressup,flag] = scioperator(b,[],'d',a.inf,a.sup);

	FlagReading(flag);
	
	if (prod(size(resinf))==0)

		c=[];

	else

		c=interval(resinf,ressup);

		if d ~= 0
			c(ii) = #(%nan,%nan);

		end

	end	

endfunction

function c=%s_q_interM(a,b)

	[a,b] = ArithReading(a,b);

	if lengthM(b) > lengthM(a) & lengthM(a) == 1 & a == 0

		c = ones(b.inf);
		c(:) = #(%nan,%nan);

		return;
	
	end

	ii = [];
	d = 0;

	for i=1:lengthM(b)
	
		if a(i) == 0
		
			ii = [i,ii];
			a(i) = 1;
			d=d+1;
		end
	end

	[resinf,ressup,flag] = scioperator(b.inf,b.sup,'d',a,[]);

	FlagReading(flag);
	
	if (prod(size(resinf))==0)

		c=[];

	else

		c=interval(resinf,ressup);

		if d ~= 0
			c(ii) = #(%nan,%nan);
		end

	end	

endfunction

// List returning for division with interval including zero

function c=extdiv(a,b)

	if ~in(0,b) then

		c = list(a/b);

	else
		
		c = list(a/#(b.inf,0),a/#(0,b.sup));

	end

endfunction

// System feedback for intervals

function c=%interM_v_interM(a,b)

	NIY();

endfunction

function c=%interM_v_s(a,b)

	NIY();

endfunction

function c=%s_v_interM(a,b)

	NIY();

endfunction

// Kroneker divisions  for intervals

function c=%interM_y_interM(a,b)

	NIY();

endfunction

function c=%interM_y_s(a,b)

	NIY();

endfunction

function c=%s_y_interM(a,b)

	NIY();

endfunction

function c=%interM_z_interM(a,b)

	NIY();

endfunction

function c=%interM_z_s(a,b)

	NIY();

endfunction

function c=%s_z_interM(a,b)

	NIY();

endfunction

// Power functions for intervals

function c=%interM_p_interM(a,b)

	NIY();

endfunction

function c=%interM_p_s(a,n)

[a,a] = ArithReading(a,a);

	if lengthM(n) == 1 then

	if n(1) == 1 then

		c = a;

	elseif n(1) == 0 then

		for i=1:lengthM(a)

		c(i) = interval(1);

		end

		c = matrix(c,size(a));

	elseif strictin(abs(n(1)),#(0,1))  then

		if lengthM(a) == 1

			if n(1) < 0

			c = 1./(a^(-round(1./n))) ;

			elseif round(1./n) == 1

			c = a;

			else

			[resinf,ressup,flag] = scipower(a.inf,a.sup,'R',round(1./n));
							
			FlagReading(flag);
	
			if (prod(size(resinf))==0)

			c=[];

			else

			c=interval(resinf,ressup);

			end	

			end
		
		else
	
			if (min(a.dims) == 1)
	
			c = a.^(n(1))

			else
		
			disp('For the moment, try with .^ !!');

			NIY();

			end 

		end

	else

		if round(n) ==  n then
	
		c = %interM_p_i(a,n);

		else
			
		NIY();

		end

	end

	else

	if lengthM(a) == 1

	for i=1:lengthM(n)
	
	c(i)=a(1)^n(i);

	end
		
	c = matrix(c,size(n));
		
	else

	NIY();
	end
			
	end

endfunction

//Integer power :

function  c=%interM_p_i(a,n)


[a,a] = ArithReading(a,a);

if lengthM(n) == 1 then	
	
if n(1) < 0

c = 1/(a^(-n(1)));

else

if (min(a.dims) ~= 1)

    [resinf,ressup,flag] = scipower(a.inf,a.sup,'P',n(1));
					
    FlagReading(flag);

    if (prod(size(resinf))==0)

       c=[];

    else

       c=interval(resinf,ressup);

    end	
else

   c=a.^n(1);
end
		
end
else

if lengthM(a)	== 1 then
			
for i=1:lengthM(n)
	
c(i)=a(1)^n(i);

end
		
c = matrix(c,size(n));

else

flag = -1;
FlagReading(flag);
end

end

endfunction

// Square root :

function c=%interM_sqrt(varargin)

	c = varargin(1).^(0.5);

endfunction

// Element by element Power functions for intervals

function c=%interM_j_interM(a,b)

	NIY();

endfunction

// Element by element root functions

function c=%interM_j_s(a,n)

[a,a] = ArithReading(a,a);

if lengthM(n) == 1 then
		
if n(1) == 1 then

c = a;

elseif n(1) == 0 then

c = interval(ones(size(a)));

elseif strictin(abs(n(1)),#(0,1))  then
 
if n(1) < 0

c = 1./(a.^(round(-1./n))) ;

else

[resinf,ressup,flag] = scipower(a.inf,a.sup,'R',round(1./n));
							
FlagReading(flag);
	
if (prod(size(resinf))==0)

c=[];

else

c=interval(resinf,ressup);

end	

end

else

if round(n) ==  n then
	
c = %interM_j_i(a,n);

else
			
NIY();

end

end

else

if size(a) == size(n) then

for i=1:lengthM(a)

c(i)=a(i)^n(i);

end

c = matrix(c,size(a));
		

elseif lengthM(a) ==1
			
for i=1:lengthM(n)
	
c(i)=a(1).^n(i);

end
	
c = matrix(c,size(n));	
	
else
	
flag = -1;

FlagReading(flag);

end
			
end

endfunction

//Element by element power functions

function  c=%interM_j_i(a,n)

[a,a] = ArithReading(a,a);

if lengthM(n) == 1 then	

   if n(1) < 0

      for i=1:lengthM(a)
	
     c(i) = 1/(a(i)^(-n));
	
     end
		
     c = matrix(c,size(a));

   else

   [resinf,ressup,flag] = scipower(a.inf,a.sup,'E',n(1));
						
   FlagReading(flag);

   if (prod(size(resinf))==0)

      c=[];

    else

      c=interval(resinf,ressup);

   end
	
end

else
	
if lengthM(a) ==1

for i=1:lengthM(n)

c(i)=a(1).^n(i);

end

c = matrix(c,size(n));	

else

flag = -1;

FlagReading(flag);

end

end

endfunction	

// Unary functions
//---------------

//Absolute value

function c=%interM_abs(a)

	[a,a] = ArithReading(a,a);

	[resinf,ressup,flag] = sciunary(a.inf,a.sup,'1');

	FlagReading(flag);

	if (prod(size(resinf))==0)

		c=[];

	else

		c=interval(resinf,ressup);

	end
		
endfunction 

//Elementwise exponential


function c=%interM_exp(a)

	[a,a] = ArithReading(a,a);

	[resinf,ressup,flag] = sciunary(a.inf,a.sup,'2');

	FlagReading(flag);

	if (prod(size(resinf))==0)

		c=[];

	else

		c=interval(resinf,ressup);

	end
		
endfunction 

//Elementwise neperian logarithm

function c=%interM_log(a)

	[a,a] = ArithReading(a,a);

	for i=1:lengthM(a)
		
		if a(i)<0
				
			disp("Invalid domain for an interval operation");
			disp("See help I4S_Arith_Func");
			error('Invalid domain for an interval operation',20012);

		end
	end

	[resinf,ressup,flag] = sciunary(a.inf,a.sup,'3');

	FlagReading(flag);

	if (prod(size(resinf))==0)

		c=[];

	else

		c=interval(resinf,ressup);

	end
		
endfunction

//Elementwise 10th based logarithm

function c=%interM_log10(a)

	[a,a] = ArithReading(a,a);

	for i=1:lengthM(a)
		
		if a(i)<0

			disp("Invalid domain for an interval operation");
			disp("See help I4S_Arith_Func");
			error('Invalid domain for an interval operation',20012);

		end
	end
	

	[resinf,ressup,flag] = sciunary(a.inf,a.sup,'4');

	FlagReading(flag);

	if (prod(size(resinf))==0)

		c=[];

	else

		c=interval(resinf,ressup);

	end
		
endfunction 

//Elementwise sinus

function c=%interM_sin(a)

	[a,a] = ArithReading(a,a);

	[resinf,ressup,flag] = sciunary(a.inf,a.sup,'a');

	FlagReading(flag);

	if (prod(size(resinf))==0)

		c=[];

	else

		c=interval(resinf,ressup);

	end
		
endfunction 


//Elementwise cosinus

function c=%interM_cos(a)

	[a,a] = ArithReading(a,a);

	[resinf,ressup,flag] = sciunary(a.inf,a.sup,'b');

	FlagReading(flag);

	if (prod(size(resinf))==0)

		c=[];

	else

		c=interval(resinf,ressup);

	end
		
endfunction 

//Elementwise tangente

function c=%interM_tan(a)

	[a,a] = ArithReading(a,a);

	for i=1:lengthM(a)

		if in(0,cos(a(i)))

			disp("Invalid domain for an interval operation");
			disp("See help I4S_Arith_Func");
			error('Invalid domain for an interval operation',20012);

		end

	end

	
	[resinf,ressup,flag] = sciunary(a.inf,a.sup,'c');

	FlagReading(flag);

	if (prod(size(resinf))==0)

		c=[];

	else

		c=interval(resinf,ressup);

	end
	
	
endfunction 

// Elementwise cotangente

function [t] = cotg(z)

	if typeof(z) == "interM"

		a = z;

		[a,a] = ArithReading(a,a);

		for i=1:lengthM(a)

			if in(0,sin(a(i)))

			disp("Invalid domain for an interval operation");
			disp("See help I4S_Arith_Func");
			error('Invalid domain for an interval operation',20012);

			end

		end

		[resinf,ressup,flag] = sciunary(a.inf,a.sup,'d');

		FlagReading(flag);

		if (prod(size(resinf))==0)

			c=[];

		else

			c=interval(resinf,ressup);

		end
		
		t = c;

	else
//Copy of internal scilab function

	if type(z)<>1 then error(53,1),end
  	t = 1 ./tan(z);

	end

endfunction

// Elementwise arcsinus

function c=%interM_asin(a)

	[a,a] = ArithReading(a,a);

	for i=1:lengthM(a)

		if ~in(a(i),#(-1,1))

			disp("Invalid domain for an interval operation");
			disp("See help I4S_Arith_Func");
			error('Invalid domain for an interval operation',20012);

		end

	end	

	[resinf,ressup,flag] = sciunary(a.inf,a.sup,'e');

	FlagReading(flag);

	if (prod(size(resinf))==0)

		c=[];

	else

		c=interval(resinf,ressup);

	end
		
endfunction

// Elementwise arcosinus

function c=%interM_acos(a)

	[a,a] = ArithReading(a,a);

	for i=1:lengthM(a)

		if ~in(a(i),#(-1,1))
			
			disp("Invalid domain for an interval operation");
			disp("See help I4S_Arith_Func");
			error('Invalid domain for an interval operation',20012);
		end

	end
	
	[resinf,ressup,flag] = sciunary(a.inf,a.sup,'f');

	FlagReading(flag);

	if (prod(size(resinf))==0)

		c=[];

	else

		c=interval(resinf,ressup);

	end
		
endfunction 


// Elementwise arctangente

function c=%interM_atan(a)

	[a,a] = ArithReading(a,a);

	[resinf,ressup,flag] = sciunary(a.inf,a.sup,'g');

	FlagReading(flag);

	if (prod(size(resinf))==0)

		c=[];

	else

		c=interval(resinf,ressup);

	end
		
endfunction 


// Elementwise hyperbolic sinus

function c=%interM_sinh(a)

	[a,a] = ArithReading(a,a);

	[resinf,ressup,flag] = sciunary(a.inf,a.sup,'i');

	FlagReading(flag);

	if (prod(size(resinf))==0)

		c=[];

	else

		c=interval(resinf,ressup);

	end
		
endfunction 


function [t] = sinh(z)

	if typeof(z) == "interM"

		t = %interM_sinh(z);

	else
// Copy of internal scilab macro

	if type(z)<>1 then error(53,1),end
  	if isreal(z) then
    	 t = imag(sin(imult(z)))
 	 else
    	 t = -imult(sin(imult(z)))
  	end
	
	end

endfunction

// Elementwise hyperbolic cosinus

function c=%interM_cosh(a)

	[a,a] = ArithReading(a,a);

	[resinf,ressup,flag] = sciunary(a.inf,a.sup,'j');

	FlagReading(flag);

	if (prod(size(resinf))==0)

		c=[];

	else

		c=interval(resinf,ressup);

	end
		
endfunction 

function [t] = cosh(z)

	if typeof(z) == "interM"

		t = %interM_cosh(z);

	else
// Copy of internal scilab macro

	if type(z)<>1 then error(53,1),end
  	if isreal(z) then
    	y = exp(abs(z)) ; t = 0.5*(y + 1 ./y)
  	else
    	t = cos(imult(z))
  	end
	
	end

endfunction

// Elementwise hyperbolic tangente

function c=%interM_tanh(a)

	[a,a] = ArithReading(a,a);

	[resinf,ressup,flag] = sciunary(a.inf,a.sup,'k');

	FlagReading(flag);

	if (prod(size(resinf))==0)

		c=[];

	else

		c=interval(resinf,ressup);

	end
		
endfunction 

function [t] = tanh(z)

	if typeof(z) == "interM"

		t = %interM_tanh(z);

	else
// Copy of internal scilab macro

	if type(z)<>1 then error(53,1),end
  	if isreal(z) then
     	t = imag(tan(imult(z)))
  	else
     	t = -imult(tan(imult(z)))
 	 end
	
	end

endfunction

// Elementwise hyperbolic cotangente

function c=%interM_coth(a)

	[a,a] = ArithReading(a,a);

	for i=1:lengthM(a)

		if in(0,a(i))

			disp("Invalid domain for an interval operation");
			disp("See help I4S_Arith_Func");
			error('Invalid domain for an interval operation',20012);

		end

	end
	[resinf,ressup,flag] = sciunary(a.inf,a.sup,'l');

	FlagReading(flag);

	if (prod(size(resinf))==0)

		c=[];

	else

		c=interval(resinf,ressup);

	end
		
endfunction 

function [t] = coth(x)

	if typeof(x) == "interM"

		t = %interM_coth(x); 

	else
// Copy of internal scilab macro
	
	if type(x)<>1 then error(53,1),end
  	t=exp(x);
  	t=(t-ones(x)./t).\(t+ones(x)./t)

	end

endfunction

// Elementwise hyperbolic arcsinus

function c=%interM_asinh(a)

	[a,a] = ArithReading(a,a);
	
	[resinf,ressup,flag] = sciunary(a.inf,a.sup,'m');

	FlagReading(flag);

	if (prod(size(resinf))==0)

		c=[];

	else

		c=interval(resinf,ressup);

	end
		
endfunction 

function [t] =  asinh(x)

	if typeof(x) == "interM"

		t = %interM_asinh(x);

	else
// Copy of internal scilab macro
	
	if type(x)<>1 then error(53,1),end
  	if isreal(x) then
    	t = imag(asin(imult(x)))
  	else
    	t = -imult(asin(imult(x)))
  	end
	
	end

endfunction

// Elementwise hyperbolic arccosinus

function c=%interM_acosh(a)

	[a,a] = ArithReading(a,a);

	for i=1:lengthM(a)

		if ~a(i)>1

			disp("Invalid domain for an interval operation");
			disp("See help I4S_Arith_Func");
			error('Invalid domain for an interval operation',20012);

		end

	end

	[resinf,ressup,flag] = sciunary(a.inf,a.sup,'n');

	FlagReading(flag);

	if (prod(size(resinf))==0)

		c=[];

	else

		c=interval(resinf,ressup);

	end
		
endfunction 

function [t] = acosh(z)

	if typeof(z) == "interM"

		t = %interM_acosh(z);

	else
// Copy of internal scilab macro
	
	if type(z)<>1 then error(53,1),end
  	if isreal(z) then
    	if min(z) < 1 then
     
      	u = acos(z)
      	t = 2*(0.5 - bool2s(imag(u)>0)).*imult(u)
   	 else
     
      	t = imag(acos(z))
    	end
 	 else
    	u = acos(z)
    	t = 2*(0.5 - bool2s(imag(u)>0)).*imult(u)
 	 end

	end

endfunction

// Elementwise hyperbolic arctangente

function c=%interM_atanh(a)

	[a,a] = ArithReading(a,a);

	for i=1:lengthM(a)

		if ~strictin(a(i),#(-1,1))

			disp("Invalid domain for an interval operation");
			disp("See help I4S_Arith_Func");
			error('Invalid domain for an interval operation',20012);

		end

	end
	[resinf,ressup,flag] = sciunary(a.inf,a.sup,'o');

	FlagReading(flag);

	if (prod(size(resinf))==0)

		c=[];

	else

		c=interval(resinf,ressup);

	end
		
endfunction 

function [t] = atanh(z)

	if typeof(z) == "interM"

		t = %interM_atanh(z);

	else
// Copy of internal scilab macro

	if type(z)<>1 then error(53,1),end

  	if isreal(z) then
    	if max(abs(z)) > 1 then 

      	t = imult(atan(-imult(z)))
   	 	else	
     
      	t= -imag(atan(-imult(z)))
    	end
  	else
   	 t = imult(atan(-imult(z)))
 	 end
	
	end

endfunction

// INTERVAL LINEAR SOLVING

function X = I4Slinearsolve(varargin)

	// degenerate intervals for real inputs 

	if (typeof(varargin(1)) == "constant")
	
		varargin(1) = #(varargin(1));

	end

	if (typeof(varargin(2)) == "constant" )
	
		varargin(2) = #(varargin(2));

	end
	
	if lengthM(varargin) >= 3
	if typeof(varargin(3)) == "constant"

		varargin(3) = #(varargin(3));

	end
	end	

	// Check input dimensions and type

	if (typeof(varargin(1)) == "interM" & typeof(varargin(2)) == "interM" )

		if size(varargin(1),'r') ~= size(varargin(1),'c')

			disp('bad dimensions ! Only square matrix for the moment');
			flag = -3;
			FlagReading(flag);
	
		end
	
		if (size(varargin(2),'r') ~= size(varargin(1),'c') | size(varargin(2),'c') ~= 1)

			disp('bad dimensions !');
			flag = -1;
			FlagReading(flag);

		end

	else
		
		disp(' Bad input type in a linear system');
		error('Bad input type in a linear system',30017);

	end
	
	if lengthM(varargin) >= 3
	if typeof(varargin(3)) == "interM" 

		if (size(varargin(3),'c') ~= size(varargin(2),'c') ..
		 | size(varargin(3),'c') ~= size(varargin(2),'c'))
			
			flag = -1;
			FlagReading(flag);

		end

	end
	end
	// Check non infinity intervals
	
	a = varargin(1);
	b = varargin(2);
	
	[a,b] = ArithReading(a,b);
	
	T1 = ( abs(a.inf) == %inf | abs(a.sup) == %inf );
	T2 = ( abs(b.inf) == %inf | abs(b.sup) == %inf );

	ii1 = find(T1);
	ii2 = find(T2);
	
	if ii1 ~= [] | ii2 ~= []
		disp ('function/operation does not handle infinite value');
		error ('function/operation does not handle infinite value', 20008);
	end

	// Default execution

	if lengthM(varargin) == 2
	
	if (typeof(varargin(1)) == "interM" & typeof(varargin(2)) == "interM")

		if (det(mid(a))) <> 0

		[resinf,ressup,flag] = sciILSR(a.inf,a.sup,b.inf,b.sup,[],[],'b');

		else

		flag = -1;
		resinf = [];
	
		end

		// Flag adaptation to FlagReading function

		if (flag == -1 )
			flag = -5 ;
		end

		FlagReading(flag);
	
		if (prod(size(resinf))==0)

			X = [];

		else

			X = interval(resinf,ressup);
		end	
		
		return;
	end
	end

	// Methods choice	
	
	if ( lengthM(varargin) == 3 & typeof(varargin(3)) == "string" )

		mc = varargin(3);

		if mc == "GE"

			[resinf,ressup,flag] = sciILSR(a.inf,a.sup,b.inf,b.sup,[],[],'e');
		
			// Finish before the end of the algorithm

			if flag == -1
				disp('Division by zero occured during algorithm !!');
				resinf = -%inf * ones(b.inf);
				ressup = %inf * ones(b.sup);
				flag = 1;
			end

		elseif mc == "B"

			[resinf,ressup,flag] = sciILSR(a.inf,a.sup,b.inf,b.sup,[],[],'f');
		
			// Finish before the end of the algorithm
		
			if flag == -1
				disp('Division by zero occured during algorithm !!');
				resinf = -%inf * ones(b.inf);
				ressup = %inf * ones(b.sup);
				flag = 1;
			end

	
		elseif mc == "LU"

			[resinf,ressup,flag] = sciILSR(a.inf,a.sup,b.inf,b.sup,[],[],'h');
		
		elseif mc == "HB"

			[resinf,ressup,flag] = sciILSR(a.inf,a.sup,b.inf,b.sup,[],[],'i');
	
			// Flag adaptation to FlagReading function
			if flag == 0
				flag = -7;
			end
		
			// Finish before the end of the algorithm

			if flag == -1
				disp('Empty solution !');
				resinf = -%nan * ones(b.inf);
				ressup = %nan * ones(b.sup);
				flag = 1;
			end
				

		elseif mc == "PGE"

			[resinf,ressup,flag] = sciILSR(a.inf,a.sup,b.inf,b.sup,[],[],'l');
		
				// Flag adaptation to FlagReading function

				if flag == -3
					flag = -6;
				end
			
				if flag == -5
					flag = -3;
				end

				// Finish before the end of the algorithm

				if flag == -1
					disp('Division by zero occured during algorithm !!');
					resinf = -%inf * ones(b.inf);
					ressup = %inf * ones(b.sup);
					flag = 1;
				end

		elseif mc == "PB"

			[resinf,ressup,flag] = sciILSR(a.inf,a.sup,b.inf,b.sup,[],[],'m');

				// Flag adaptation to FlagReading function
				if flag == -3
					flag = -6;
				end
		
				if flag == -5
					flag = -3;
				end

				// Finish before the end of the algorithm
				if flag == -1
					disp('Division by zero occured during algorithm !!');
					resinf = -%inf * ones(b.inf);
					ressup = %inf * ones(b.sup);
					flag = 1;
				end

		elseif mc == "PLU"

			[resinf,ressup,flag] = sciILSR(a.inf,a.sup,b.inf,b.sup,[],[],'o');
		
		elseif mc == "PHB"

			[resinf,ressup,flag] = sciILSR(a.inf,a.sup,b.inf,b.sup,[],[],'p');

				// Flag adaptation to FlagReading function
				
				if flag == -3
					flag = -6;
				end

				if flag == -5
					flag = -3;
				end

				if flag == 0
					flag = -7;
				end
		
				// Finish before the end of the algorithm

				if flag == -1
					disp('Empty solution !');
					resinf = -%nan * ones(b.inf);
					ressup = %nan * ones(b.sup);

					flag = 1;
				end

		else
			disp('Bad method choice for a linear system');
			error('Bad method choice for a linear system', 30018);

		end
			
	elseif ( lengthM(varargin) >= 4 & typeof(varargin(4)) == "string" ..
		& typeof(varargin(3)) == "interM" )

        I4Sconst = tlist(["I4Sconst","Max_Loop","Fac_improve"], 10, 0.9);

		mc = varargin(4);
	
		Xin = varargin(3);

		MF = I4Sget2C();

		dd = 0;
		
		if lengthM(varargin) == 6
			
			if varargin(6) ~= []
			if (MF(3) ~= varargin(6) | varargin(6) ~= I4Sconst(3))
			MF(3) = varargin(6);
			dd = dd+1;
			end
			else
			MF(3) = I4Sconst(3);
			dd = dd+1;
			end

			if  varargin(5) ~= []
			if (MF(2) ~= varargin(5) | varargin(5) ~= I4Sconst(2))
			MF(2) = varargin(5);
			dd = dd+1;
			end
			else
			MF(2) = I4Sconst(2);
			dd = dd+1;
			end

		end

		if (lengthM(varargin) == 5  & varargin(5) ~= [])
			
			if (MF(2) ~= varargin(5) | varargin(5) ~= I4Sconst(2))
			MF(2) = varargin(5);
			dd = dd+1;
			end
		end

		if dd ~= 0
			I4Ssend2C(MF);
		else
			I4Ssend2C(I4Sconst);
		end

		[Xin,Xin] = ArithReading(Xin,Xin);
		
		if mc == "K"

			[resinf,ressup,flag] = sciILSR(a.inf,a.sup,b.inf,b.sup,Xin.inf,Xin.sup,'c');
		
			// Finish before the end of the algorithm
			
			if flag == 0
				disp(' non intersection during iterative Krawczyk algorithm');	
				flag = 1;
			end

			if flag == -1	
				disp(' Krawczyck''s test failure');
				resinf = Xin.inf;
				ressup = Xin.sup;
				flag = 1;
			end	

		elseif mc == "KGS"

			[resinf,ressup,flag] = sciILSR(a.inf,a.sup,b.inf,b.sup,Xin.inf,Xin.sup,'d');
		
			// Finish before the end of the algorithm

			if flag == 0
				disp(' non intersection during iterative Krawczyk algorithm');	
				flag = 1;
			end

			if flag == -1	
				disp(' Krawczyck''s test failure');
				resinf = Xin.inf;
				ressup = Xin.sup;
				flag = 1;
			end	
		
		elseif mc == "GS"

			[resinf,ressup,flag] = sciILSR(a.inf,a.sup,b.inf,b.sup,Xin.inf,Xin.sup,'g');
		
			// Finish before the end of the algorithm

			if flag == -1
				disp(' no more intersection during iterative Gauss Seidel algorithm');	
				flag = 1;
			end

		elseif mc == "PK"

			[resinf,ressup,flag] = sciILSR(a.inf,a.sup,b.inf,b.sup,Xin.inf,Xin.sup,'j');

				// Flag adaptation to FlagReading function
			
				if flag == -3
					flag = -6;
				end
		
				if flag == -5
					flag = -3;
				end

				// Finish before the end of the algorithm
			
				if flag == 0
					disp(' non intersection during iterative Krawczyk algorithm');	
					flag = 1;
				end

				if flag == -1	
					disp(' Krawczyck''s test failure');
					resinf = Xin.inf;
					ressup = Xin.sup;
					flag = 1;
				end	

		elseif mc == "PKGS"

			[resinf,ressup,flag] = sciILSR(a.inf,a.sup,b.inf,b.sup,Xin.inf,Xin.sup,'k');

				// Flag adaptation to FlagReading function
			
				if flag == -3
					flag = -6;
				end

				if flag == -5
					flag = -3;
				end

				// Finish before the end of the algorithm
			
				if flag == 0
					disp(' non intersection during iterative Krawczyk algorithm');	
					flag = 1;
				end

				if flag == -1	
					disp(' Krawczyck''s test failure');
					resinf = Xin.inf;
					ressup = Xin.sup;
					flag = 1;
				end	
		
		elseif mc == "PGS"

			[resinf,ressup,flag] = sciILSR(a.inf,a.sup,b.inf,b.sup,Xin.inf,Xin.sup,'n');

				// Flag adaptation to FlagReading function

				if flag == -3
					flag = -6;
				end

				if flag == -5
					flag = -3;
				end

				// Finish before the end of the algorithm

				if flag == -1
					disp(' no more intersection during iterative Gauss Seidel algorithm');	
					flag = 1;
				end

		else
			disp('Bad method choice for a linear system');
			error('Bad method choice for a linear system', 30018);

		end	
	else

		disp(' Bad input type in a linear system');
		error('Bad input type in a linear system',30017);

	end		
	
	// Flag adaptation to FlagReading function

	if  flag == -2 
		flag = -6;
	end

	if flag == -4
		flag = -5;
	end

	FlagReading(flag);
	
	if (prod(size(resinf))==0)

		X = [];

	else

		X = interval(resinf,ressup);
	end
	
endfunction

function a = I4Sget2C()

    [c,d]= sciI4Svarget('0');
    b = tlist(["I4Sconst","Max_Loop","Fac_improve"]);
    b.Max_Loop = d;
    b.Fac_improve = c;
    a = b;
	  //disp(" Int4Sci constants are now :");	
	  //disp(a);
 
endfunction

function I4Ssend2C(a)

    if typeof(a) == "I4Sconst";
    if round(a.Max_Loop) ~= a.Max_Loop
      disp("Max Loop must be an integer !!");
      abort;
    end
    sciI4Svarsend('0',a.Fac_improve,a.Max_Loop);
    end
 
endfunction

function %I4Sconst_p(x)
	
	st(1) = "---------------------------------";		
	st(2) = "---------------------------------";
	st(3) = "Max_Loop : " + string(x.Max_Loop);
  st(4) = st(2);
  st(5) = "Fac_improve : " + string(x.Fac_improve);
  st(6) = st(2);
	st= st(:);
	disp(st);

endfunction

