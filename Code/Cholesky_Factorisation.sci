/*
Etideur:
    Jinshan GUO
Objecitf:
    Fonction à réaliser l'algorithme de Cholesky pour le calcul direct en factorisation
Principe
    A = C * C^T 
Containtes:
    A est une matrice symétrique définie positive
    C est une matrice triangulaire inférieure et C(i,i)>0
    C^T est matrice transposée de C, triangulaire supérieure
Valeur retour:
    C
*/
function [C,CT]=Cholesky_Factorisation(A)
    n = size(A)(1); //Nombre de ligne de matrice
    C = zeros(n,n);
    for j=1:n
        s = A(j,j) - C(j,1:j) * C(j,1:j)';
        if abs(s) < %eps then
            error("A non SDP");
        end
        C(j,j) = sqrt(s);
        /*for i=j+1:n
            C(i,j) = 1/C(j,j) * (A(i,j) - C(i,1:j) * C(j,1:j)');
        end
        */
        C(j+1:n,j) = 1/C(j,j) * (A(j+1:n,j) - C(j+1:n,1:j) * C(j,1:j)');
    end
    CT = C';
endfunction

