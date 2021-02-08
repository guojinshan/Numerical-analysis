/*
Créateur: Jinshan GUO et Anais Debureaux
*/

//===========Exo1-solsup===========
function [x]=solsup(U,b)
    n = size(U)(1); //Nombre lignes de matrice
    x = zeros(n,1);
    x(n) = 1/U(n,n) * b(n);
    for i=n-1:-1:1
        if abs(U(i,i)) < %eps then
            error("Matrice U est non inversible");
        end
        x(i) = 1/U(i,i) * (b(i) - U(i,i+1:n) * x(i+1:n));
    end
endfunction

//===========Exo2-trisup===========
function [U,e]=trisup(A,b)
    n = size(A)(1); //Nombre de ligne de matrice
    for k=1:n-1
        if abs(A(k,k)) < %eps then
            error('Matrice est non inversible pendant élimination.');
        else
            for i = (k+1):n
                coef = A(i,k)/A(k,k);
                b(i) = b(i) - coef*b(k);
                A(i,k) = 0;
                A(i,k+1:n) = A(i,k+1:n) - coef*A(k,k+1:n);
            end
        end      
    end
    if abs(A(n,n)) < %eps then
        error("Le dernier pivot est null, Matrice non inversible pendant l''élimination.");
    end
    U = A;
    e = b;
endfunction

//===========Exo3-resolG===========
function [x]=resolG(A,b)
    [U,e]=trisup(A,b);
    [x]=solsup(U,e);
endfunction

//===========Exo4===========
/*
Question a  
 F_i(ligne i de F) = E_i(ligne i de E) - C_i(ligne i de C) * V
 F_j(colonne j de F) = E_j(colonne j de E) - C * V_i(colonne i de V) 
 avec i et j allant de i à p
*/

//Question b
function [F]=resol(E,c,v)
    [U,e]=trisup(A,b);
    [x]=solsup(U,e);
endfunction

//===========Exo5-LU===========
function [L,U] = LU(A)
    n = size(A)(1); //Nombre de ligne de matrice
    L = eye(n,n);
    for k=1:n-1
        if abs(A(k,k)) < %eps then
            error('Matrice est non inversible pendant élimination.');
        else
            for i = (k+1):n
                coef = A(i,k)/A(k,k);
                L(i,k) = coef; 
                A(i,k) = 0;
                A(i,k+1:n) = A(i,(k+1):n) - coef*A(k,(k+1):n);
            end
        end      
    end
    if abs(A(n,n)) < %eps then
        error("Le dernier pivot est null, Matrice non inversible pendant l''élimination.");
    end
    U = A;
endfunction

//===========Exo6-solinf===========
function [x]=solinf(L,b)
    n = size(L)(1); //Nombre lignes de matrice
    x = zeros(n,1);
    x(1) = b(1)/L(1,1);
    for i=2:n 
        if abs(L(i,i)) < %eps then
            error("Matrice L est non inversible");
        end
        x(i) = 1/L(i,i)*(b(i)-L(i,1:i-1)*x(1:i-1));
    end
endfunction

//===========Exo7-resolLU===========
function [x] = resolLU(A,b)
    //Calculer L et U
    [L,U] = LU(A) 
    //Calculer y
    [y] = solinf(L,b);
    //Calculer x
    [x]=solsup(U,y);
endfunction

//===========Exo8-inverse===========
function [B]=inverse(A)
    n = size(A)(1); //Nombre lignes de matrice
    I = eye(n,n);
    B = [];
    B = [B,resolLU(A,I(:,1:n))];
endfunction





