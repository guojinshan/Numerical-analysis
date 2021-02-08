function y= LagrangePol(x, xj, i)
       // Computes L_i(x) with interpolation points xj 
       // at points x
       // Entry : x=vector of points where to evaluate L_i
       numer = ones(x); // init numer with ones of size of x
       denom = 1.0;
       l = length(xj);
       setindex = setdiff(1:l, i); // extract i from the index
       disp(setindex)
       //
       for j=setindex // loop j<>i
           numer = numer .* (x-xj(j));
           denom = denom .* (xj(i)-xj(j));
       end;
       y = numer / denom;
    endfunction
