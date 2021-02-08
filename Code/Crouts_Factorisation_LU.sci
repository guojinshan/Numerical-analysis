/*
Etideur:
    Jinshan GUO
Objecitf:
    Fonction à réaliser l'algorithme de Crout pour le calcul direct en factorisation LU
Principe
    A = LU => Ax = b => LUx = b => Ly = b,Ux = y 
Containtes:
    A est une matrice carré inversible
    L est une matrice triangulaire inférieure 
    U est une matrice triangulaire supérieure dont les termes diagonaux sont 1
    b est un vecteur colonne
    y est un vecteur colonne intermédiaire
    x est un vecteur colonne objetif à obtenir
Valeur retour:
    L,U,x
*/
 
 function [L,U,x] = Crout_Factorisation_LU(A,b)
     n = size(A)(1) //Nombre de ligne de matrice
     //Initialisation de L et U
     L = zeros(n,n);
     U = eye(n,n);
     for k=1:n
         //disp(L);
         //disp(U);
         L(k:n,k) = A(k:n,k) - L(k:n,1:k) * U(1:k,k); 
         if abs(L(k,k)) < %eps then
            error('Pivot nul');
         end
         U(k,k+1:n) = 1/L(k,k) * (A(k,k+1:n) - L(k,1:k) * U(1:k,k+1:n));
     end
    //Calculer y
    exec(fullpath(pwd() + '\Système Triangulaire Inférieur_SL.sci'),-1);
    [y] = sloinf(L,b);
    //Calculer x
    exec(fullpath(pwd() + '\Système Triangulaire Supérieur_SL.sci'),-1);
    [x] = slosup(U,y);
 endfunction

