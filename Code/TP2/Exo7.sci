function [x] = resolLU(A,b)
    //Calculer L et U
    exec(fullpath(pwd() + '\Exo5.sci'),-1);
    [L,U] = LU(A)
    //Calculer y
    exec(fullpath(pwd() + '\Exo6.sci'),-1);
    [y] = solinf(L,b);
    //Calculer x
      exec(fullpath(pwd() + '\Exo1.sci'),-1);
    [x]=solsup(U,y);
endfunction
