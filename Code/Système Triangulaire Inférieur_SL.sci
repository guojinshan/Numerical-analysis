/*
Etideur:
    Jinshan GUO
Objectif:
    Fonction à resoudre le problème de Ly=b
Containtes:
    L est une matrice carrée triangulaire inférieur
    b est un vecteur colonne
Valeur retour:
    y 
*/

function [y]=sloinf(L,b)
    nb_lignes = size(L)(1); //Nombre lignes de matrice
    for i=1:nb_lignes
        if abs(L(i,i)) < %eps then
            error("Matrice L est non inversible");
        end
        sum_temp = 0;
        for j=1:i-1
            sum_temp = sum_temp + L(i,j)*y(j);
        end
        y(i) = 1/L(i,i)*(b(i)-sum_temp);
    end
endfunction
