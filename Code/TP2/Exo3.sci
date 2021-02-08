function [x]=resolG(A,b)
    exec(fullpath(pwd() + '\Exo2.sci'),-1);
    [U,e]=trisup(A,b)
    exec(fullpath(pwd() + '\Exo1.sci'),-1);
    [x]=solsup(U,e);
endfunction
