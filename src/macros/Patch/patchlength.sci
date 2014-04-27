

function R=length(varargin)
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
