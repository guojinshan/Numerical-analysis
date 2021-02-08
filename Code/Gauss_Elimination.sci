/*
Etideur:
    Jinshan GUO
Objecitf:
    Fonction à réaliser l'algorithme d'élimination de Gausse
Containtes:
    A est une matrice carré inversible
    A_hat est matrice triangulaire supérieur obtenue après l'élimination
    b est un vecteur colonne
    b_hat est vecteur colonne obtenu après l'élimination
    x est un vecteur colonne objetif à obtenir
Valeur retour:
    A_hat,b_hat,x
*/

function [A_hat,b_hat,x]=Gauss_Elimination(A,b)
    n = size(A)(1); //Nombre de ligne de matrice
    for k=1:n-1
        //disp(A); 
        //disp(b);
        if abs(A(k,k)) < %eps then
            error('Matrice est non inversible pendant élimination.');
        else
            for i = (k+1):n
                coef = A(i,k)/A(k,k);
                b(i) = b(i) - coef*b(k);
                A(i,k) = 0;
                for j=(k+1):n
                    A(i,j) = A(i,j) - coef*A(k,j);
                end
            end
        end      
    end
    if abs(A(n,n)) < %eps then
        error("Le dernier pivot est null, Matrice non inversible pendant l''élimination.");
    end
    A_hat = A;
    b_hat = b;
    PATH_STS = fullpath(pwd() + '\Système Triangulaire Supérieur_SL.sci')
    exec(PATH_STS,-1);
    [x] = slosup(A_hat,b_hat);
endfunction
