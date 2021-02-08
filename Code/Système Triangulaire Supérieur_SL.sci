/*
Etideur:
    Jinshan GUO
Objectif:
    Fonction à resoudre le problème de Ux=y
Containtes:
    U est une matrice carrée triangulaire supérieur
    y est un vecteur colonne
Valeur retour:
    x 
*/

function [x]=slosup(U,y)
    nb_lignes = size(U)(1); //Nombre lignes de matrice
    for i=nb_lignes:-1:1
        if abs(U(i,i)) < %eps then
            error("Matrice U est non inversible");
        end
        sum_temp = 0;
        for j=i+1:nb_lignes
            sum_temp = sum_temp + U(i,j)*x(j);
        end
        x(i) = 1/U(i,i)*(y(i)-sum_temp);
    end
endfunction
