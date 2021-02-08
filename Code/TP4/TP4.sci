/*
Créateur: Jinshan GUO et Anais Debureaux
*/
exec(fullpath(pwd() + '\TP4.sci'),-1);

//===========Exo1-Factorisation de Cholesky===========
function [C]=cholesky(A)
    n = size(A,1); //Nombre de ligne de matrice
    C = zeros(n,n);
    for j=1:n
        s = A(j,j) - C(j,1:j) * C(j,1:j)';
        if abs(s) < %eps then
            error("A non SDP");
        end
        C(j,j) = sqrt(s);
        C(j+1:n,j) = 1/C(j,j) * (A(j+1:n,j) - C(j+1:n,1:j) * C(j,1:j)');
    end
endfunction

function [x]=solsup(U,b)
    n = size(U,1); //Nombre lignes de matrice
    x = zeros(n,1);
    x(n) = 1/U(n,n) * b(n);
    for i=n-1:-1:1
        if abs(U(i,i)) < %eps then
            error("Matrice U est non inversible");
        end
        x(i) = 1/U(i,i) * (b(i) - U(i,i+1:n) * x(i+1:n));
    end
endfunction

function [x]=solinf(L,b)
    n = size(L,1); //Nombre lignes de matrice
    x = zeros(n,1);
    x(1) = b(1)/L(1,1);
    for i=2:n 
        if abs(L(i,i)) < %eps then
            error("Matrice L est non inversible");
        end
        x(i) = 1/L(i,i)*(b(i)-L(i,1:i-1)*x(1:i-1));
    end
endfunction

function [x]=resolchol(A,b)
    //Calcule C
    [C] = cholesky(A);
    //Calcule y
    [y] = solinf(C,b);
    //Calcule x
    [x] =  solsup(C',y);
endfunction

//===========Exo2-Matrice Hibert===========
function [H]=hilbermat(n)
    //Initialiser H
    H = zeros(n,n);
    for i=1:n
        for j=1:n
            H(i,j) = 1/(i+j-1);
        end
    end
endfunction

function condit = cond_norm_2(H)
    //Obtenir toutes les valeurs propres de H
    vp = spec(H);
    //Car H est symétrique, on simplifie le calcul
    condit = max(abs(vp)) / min(abs(vp));
endfunction

//===========Exo3-Point fixe===========
function [C]=pointfixe(g,Kmax,tol,x0)
    x(1) = x0;
    for k=1:Kmax
        x(k+1) = g(x(k));
        if norm(x(k+1) - x(k)) / norm(x(k)) < tol then
             C = x(k+1);
             return;
        else
            x(k) = x(k+1);
        end
    end
    error("Non convergece");
endfunction

//===========Exo4-Matrice creuse===========
function A2=Constuire_A2(n)
    A2 = (diag(ones(n,1)+1) + diag(ones(n-1,1)-2,-1) +diag(ones(n-1,1)-2,1)) / n;
endfunction

function A3=Constuire_A3(n)
    A3 = diag(n:-1:1);
    for i=1:n-1
       A3 = A3 + diag([1;zeros(n-1-i,1)],i) + diag([1;zeros(n-1-i,1)],-i);
    end
endfunction

