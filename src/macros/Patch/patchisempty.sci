function r =isempty(a)

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

	r = oldisempty(a);

end

endfunction