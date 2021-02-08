/*
Etideur:
    Jinshan GUO
Objecitf:
    Fonction à réaliser l'algorithme de Doolittle pour le calcul direct en factorisation LU
Principe
    A = LU => Ax = b => LUx = b => Ly = b,Ux = y 
Containtes:
    A est une matrice carré inversible
    L est une matrice triangulaire inférieure dont les termes diagonaux sont 1
    U est une matrice triangulaire supérieure
    b est un vecteur colonne
    y est un vecteur colonne intermédiaire
    x est un vecteur colonne objetif à obtenir
Valeur retour:
    L,U,x
*/
 
 function [L,U,x] = Doolittle_Factorisation_LU(A,b)
     n = size(A)(1) //Nombre de ligne de matrice
     //Initialisation de L et U
     L = eye(n,n);
     U = zeros(n,n);
     for k=1:n
         //disp(L);
         //disp(U);
         U(k,k:n) = A(k,k:n) - L(k,1:k) * U(1:k,k:n);
         if abs(U(k,k)) < %eps then
            error('Pivot nul');
         end
         L(k+1:n,k) = 1/U(k,k) * (A(k+1:n,k) - L(k+1:n,1:k) * U(1:k,k));
     end
    //Calculer y
    exec(fullpath(pwd() + '\Système Triangulaire Inférieur_SL.sci'),-1);
    [y] = sloinf(L,b);
    //Calculer x
    exec(fullpath(pwd() + '\Système Triangulaire Supérieur_SL.sci'),-1);
    [x] = slosup(U,y);
 endfunction
