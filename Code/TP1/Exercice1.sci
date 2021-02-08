function [x]=suite(N,a,b)
    x(1) = a;
    x(2) = b;
    for n=1:N-2
        x(n+2) = 9/4*x(n+1) - 1/2*x(n); 
    end
endfunction
