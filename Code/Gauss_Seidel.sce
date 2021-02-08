/*
Etideur:
    Jinshan GUO
Objecitf:
    Fonction à réaliser l'algorithme de la méthode de Gauss_Seidel pour résoudre le système linéaire Ax = b
Principe:
    A = M - N = D - E - F avec M = D + E , N = F
    M*x^(k+1) = N*x^(k) + b => x^(k+1) = M^(-1)*N*x^(k) + M^(-1)*b
Containtes:
    A est une matrice inversible
    D est une matrice de la  partie diagonale de A
    E est une matrice de la  partie triangulaire inférieure de A, avec E = -A
    F est une matrice de la  partie triangulaire supérieure de A, avec F = -A 
Valeur retour:
    x est la résolution du système
    k est nombre d'intération à converger
*/
function [x,k]=Gauss_Seidel(A, b, x0, Kmax, tol)
    /*
        x0: valeur initiale
        Kmax: nombre maximum de l'itération
        tol: tolérance de l'erreur
    */
    n = size(A,1); //Nombre de ligne de la matrice A
    if length(x0) <> n then
        error("Taille du vecteur initial incorrecte");
    end
    x = zeros(n,1);
    for k=1:Kmax
        for i=1:n
            if abs(A(i,i)) < tol then
                error("Matrice A non invesible")
            end
            s_x_k = 0
            for j = i+1:n
                s_x_k = s_x_k + A(i,j) * x0(j);
            end
            s_x_k_1 = 0
            for j = 1:i-1
                s_x_k_1 = s_x_k_1 + A(i,j) * x(j);
            end
            x(i) = (1 / A(i,i)) * (b(i) -s_x_k_1 - s_x_k ); // Calculer x^(k+1) 
        end
        if norm(x-x0)/norm(x) < tol then // Tester la convergence
            return;
        else
            x0 = x;
        end
    end
    if k == Kmax then
        error("La méthode Jacobi non convergente");
    end
endfunction

function [x,k]=Gauss_Seidel_Mat(A, b, x0, Kmax, tol)
    //Vérification: aucun terme de la diagonal de A n'est nul
    if ~and(diag(A)) then
        error("erreur: présence d''un zéro sur la diagonale de A");
    end
    // Décomposition de A = D - E- F
    D = diag(diag(A));
    E = - (triu(A) - D);
    F = - (tril(A) - D);
    // Initialisation
    x = x0;
    // Boucle itérative de résolution 
    for k = 1:Kmax
        x = inv(D - E) * (F * x + b);
        if max(abs(A*x - b)) < tol then
            return;
        end
    end
    if k == Kmax then
        error("La méthode Gauss-Seidel non convergente");
    end
endfunction

A = [6 2 3; -1 7 5; 3 -2 6]; //Matrice diagonale strictement dominante
b = [9;8;7];
x0 = [1;1;2];
Kmax = 1000;
tol = 0.0000001;

[x1,k1]=Gauss_Seidel(A, b, x0, Kmax, tol)
disp(x1);
disp(k1);

[x2,k2]=Gauss_Seidel_Mat(A, b, x0, Kmax, tol)
disp(x2);
disp(k2);
