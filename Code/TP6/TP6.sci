/*
Créateur: Jinshan GUO et Anais Debureaux
*/

//===========Exo1===========
function [A] = construct(t,tau)
    n = length(t);
    m = length(tau);
    if m<n then
        error("m inférieur à n");
    end
    //Initialisation de A
    A = zeros(m,n);
    for i=1:m
        // Obtenir position indice j
        j = place(t,tau(i))
        A(i,j) = (tau(i) - t(j+1))/(t(j) - t(j+1))
        A(i,j+1) = (tau(i) - t(j))/(t(j+1) - t(j))
    end
endfunction

//===========Exo2===========
function [x] = mcnorm(A, y)
    //Cas tau veteur ligne
    if size(y,1) < size(y,2) then
        y = y';
    end
    // Cas A carrée
    if size(A,1) == size(A,2) then
        if det(A) == 0 then
            error('A non inversible');
        else
            x = A \ y;
            //x=resolchol(A,y)
        end
    else
    // Cas A non carrée
        ATA = A' * A;
        ATy = A' * y;
        if det(ATA) == 0 then
            error('ATA non inversible');
        else
            x = ATA \ ATy
            //x=resolchol(ATA,ATy)
        end
    end
endfunction

//===========Exo3===========
function [x, k] = newton(foncjac, tol, N, x0)
    grad = 1
    k = 0
    x = x0
    [F_x0, J_x0] = foncjac(x0);
    while grad > tol & k < N
        [F, J] = foncjac(x);
        step = J \ F;
        grad = norm(F)/norm(F_x0);
        x = x - step;
        k = k + 1;
    end
    if k==N then
        error("La méthode de Newton n''a pas convergé");
    end
endfunction

function courbe(N)
    //scf();
    [x, k] = newton(foncjac, tol, N, x0);
    y_1 = x(1)^2 + x(2)^2 -4;
    y_2 = x(2) - exp(x(1));
endfunction

function [F, J] = foncjac(x)
    F = [x(1)^2 + x(2)^2 - 4; x(2) - exp(x(1))];
    J = [2*x(1) 2*x(2); -exp(x(1)) 1];
endfunction

//===========Exo4===========
function [x]=mcQR(A,y)
    [Q,R]=qr(A);
    c=Q'*y;
    n = size(A,2);
    Rt=R(1:n,:);
    c1=c(1:n);
    x=Rt\c1; //ou x=solsup (Rt,c1)
endfunction
//===========Fonction du TP précedent===========
function [j] = place(t,tau)
    n = length(t);
    if tau<t(1)| tau>t(n) then
        error("La valeur de tau est hors de portée");
    end
    if tau == t(n) then
        j = n-1; 
        return;
    end
    i_min = 1;
    i_max = n;
    while i_max - i_min > 1
        i_mil = floor((i_max + i_min) / 2);
        if tau >= t(i_mil) then
           i_min = i_mil;
        else
            i_max = i_mil;
        end
    end
    j = i_min;
endfunction

function [C] = cholesky(A)
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

function [x] = solsup(U,b)
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

function [x] = solinf(L,b)
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










