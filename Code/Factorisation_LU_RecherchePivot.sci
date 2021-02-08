/*
Etideur:
    Jinshan GUO
Objecitf:
    Fonction à réaliser l'algorithme de factorisation LU avec recherche de pivots(PA=LU) 
Containtes:
    A est une matrice carré inversible
    L est une matrice triangulaire inférieure obtenue après la factorisation
    U est une matrice triangulaire supérieure obtenue après la factorisation
    P est une matrice de permutation(avec que des 0, un seul 1 par ligne et par colonne)
    b est un vecteur colonne
    y est un vecteur colonne intermédiaire
    x est un vecteur colonne objetif à obtenir
Valeur retour:
    L,U,P,x
*/

function [L,U,P,x]=Factorisation_LU_RecherchePivot(A,b)
    n = size(A)(1) //Nombre de ligne de matrice
    //Initialisation de L et P
    L = eye(n,n);
    P = eye(n,n);
    for k=1:n-1
        //Trouver le pivot maximal
        pivot_max = A(k,k);
        nombre_ligne = k;
        for i=k+1:n
            if A(i,k) > pivot_max then
                pivot_max = A(i,k);
                nombre_ligne = i;
            end
        end
        if abs(A(nombre_ligne,k)) < %eps then
            error("Pivots sont tous nulls après sa recherche, matrice inversible");
        elseif nombre_ligne ~= k then 
            //Appliquer la permutation de lignes (k,l) à P,A,L dans les colonnes 1 à k-1
             //Permutation P
            temp_P = P(k,1:k-1);
            P(k,1:k-1) = P(nombre_ligne,1:k-1);
            P(nombre_ligne,1:k-1) = temp_P;
             //Permutation A
            temp_A = A(k,1:k-1);
            A(k,1:k-1) = A(nombre_ligne,1:k-1);
            A(nombre_ligne,1:k-1) = temp_A;
             //Permutation L
            temp_L = L(k,1:k-1);
            L(k,1:k-1) = L(nombre_ligne,1:k-1);
            L(nombre_ligne,1:k-1) = temp_L;
        end
        //Effectation d'éliminqtion de Gauss sur A permutée
        for i=k+1:n
            L(i,k) = A(i,k)/A(k,k);
            for j=k:n
                A(i,j) = A(i,j) - L(i,k) * A(k,j);
            end
        end
    end
    if abs(A(n,n))< %eps then
        error("Dernier pivot est null");
    end
    U = A;
    //Calculer y
    exec(fullpath(pwd() + '\Système Triangulaire Inférieur_SL.sci'),-1);
    [y] = sloinf(L,b);
    //Calculer x
    exec(fullpath(pwd() + '\Système Triangulaire Supérieur_SL.sci'),-1);
    [x] = slosup(U,y);
endfunction
