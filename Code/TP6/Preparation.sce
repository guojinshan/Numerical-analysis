/*
Créateur: Jinshan GUO et Anais Debureaux
*/

/*
On veut minimiser E(x1, x2, x3) = \sum_{i}^{m}(p(tau_{i})-y_{i})^2
Cela revient à trouver le minimum de \norm{Ax-y}_{2}^{2}, avec :
 A = [1 t1 t1^2; 1 t2 t2^2; ...; 1 tm tm^2]
 x = [x0; x1; x2]
 y = [y1; y2; y3]
 n = 3 = nombre d'inconnus : x1 x2 x3
 */

disp("Resulat_Exo1(a)_n=", 3);

function [A] = constrpoly(tau)
    //Vérifier la taille de tau
    if length(tau) < 3 then
        error("Taille de tau inférieure à 3");
    end
    //Cas tau vecteur ligne (nombre de lignes < nombre de colonnes)
    if size(tau,1) < size(tau,2) then
        tau = tau'; //on transpose le vecteur tau
    end
    // Construire la matrice de Vandermonde
    A = [tau.^0 tau.^1 tau.^2];
endfunction

function [x] = mcnorm(A, y)
    //Cas tau vecteur ligne
    if size(y,1) < size(y,2) then
        y = y';
    end
    // Cas A carrée
    if size(A,1) == size(A,2) then
        if det(A) == 0 then
            error('A carrée et non inversible');
        else
            disp ("A carrée est inversible");
            x = A \ y;
        end
    else
    // Cas A non carrée
        disp ("x est une solution au sens des moindres carrés : norm(A*x-b) est minimale");
        ATA = A' * A;
        ATy = A' * y;
        x = ATA \ ATy;
        /*
        if rank(A) == size(A,2) then
            //A est de rang maximal (colonnes linéairement indépendantes)
            //disp ("Solution unique "); 
            x = ATA \ ATy;
        else
            //disp ("Solution pas unique ");
            x = ATA \ ATy;
        end
        */
    end
endfunction

tau = [0 1 2];
y = [1 3 7];
A = constrpoly(tau);
x = mcnorm(A, y);
disp("Resulat_Exo2(b)_A=", A);
disp("Resulat_Exo2(b)_x=", x);
p = A * x;
disp("Resulat_Exo2(b)_p=", p);
disp("Le résultat semble correcte car A est inversible et on vérifie : la somme d''une ligne de A = la même ligne de p");

t = linspace(-1,7,30);
A_t = constrpoly(t);
p_t = A_t * x;
plot(t, p_t, 'r'); //parabole
scatter(tau, y, marker=4); //points de coordonnees(tau(i); yi)
xgrid; //ajoute une grille 
xlabel('t','fontsize',3);
ylabel('p(t)','fontsize',3);
title("Application1",'fontsize',4);
legend('Parabole obtenue','Point coordonnées',1);

tau = linspace(0, 2, 9);
y = [1 1.7 1.95 1.8 3. 3.6 4.45 5.9 6.6];
A = constrpoly(tau);
disp(A);
x = mcnorm(A, y);
disp("Resulat_Exo2(c)_x=", x);
scf(); //crée une nouvelle fenêtre pour le graphique
p = A * x;
plot(tau, p, 'r');
scatter(tau, y, marker=4);
xgrid;
xlabel('tau','fontsize',3);
ylabel('p(tau)','fontsize',3);
title("Application2", 'fontsize',4);
legend('Parabole obtenue','Point coordonnées',1);


