/*
Créateur: Jinshan GUO et Anais Debureaux
*/

exec(fullpath(pwd() + '\Routine.sci'),-1);
//===========Exo1-solsup===========
U = [1 2 3;0 4 8;0 0 5];
b = [6 16 15];
[x]=solsup(U,b);
disp("Exo1-x=",x);
//Vérification par Scilab
//disp("Exo1-x_resultat_scilab",U\b);

//===========Exo2-trisup===========
A = [3 1 2;3 2 6;6 1 -1];
b = [2;1;4];
[U,e]=trisup(A,b)
disp("Exo2-U=",U);
disp("Exo2-e=",e);

//===========Exo3-resolG===========
A1 = [1 2 3;5 2 1;3 -1 1];
b1 = [5;5;6];
[x1]=resolG(A1,b1);
disp("Exo3-x1=",x1);
disp("Exo3-resultat_scilab_x1=",A1\b1);

A2 = [2 1 5;1 2 4;3 4 10];
b2 = [5;5;6];
//L'éxecution s'arrête quand l'erreur apparaît car A n'est pas inversible sur le dernier pivot, on ajoute un commentaire pour assurer que le script fonctionne sans bloqué.
disp("Exo3-x2=","Bloqué car le dernier pivot est null, matrice non inversible");
//[x2]=resolG(A2,b2); 
//disp("Exo3-x2=",x2); 
disp("Exo3-x2_resultat_scilab=",A2\b2);

//===========Exo4===========



//===========Exo5-LU===========
A = [3 1 2;3 2 6;6 1 -1];
[L,U] = LU(A);
disp("Exo5-L=",L);
disp("Exo5-U=",U);

//===========Exo6-solinf===========
L = [1 0 0;2 3 0;1 4 -1];
b = [1;8;10];
[x]=solinf(L,b);
disp("Exo6-x=",x);
//Vérification par Scilab
//disp("Exo6-x_resultat_scilab",L\b);

//===========Exo7-resolLU===========
A = [1 2 3;5 2 1;3 -1 1];
[x] = resolLU(A,b);
disp("Exo7-x=",x);
//Vérification par Scilab
//disp("Exo7-x_resultat_scilab",A\b);

//===========Exo8-inverse===========
A = [1 2 3;5 2 1;3 -1 1];
[B]=inverse(A);
disp("Exo8-A_Inverse=",B);
//Vérification par Scilab
disp("Exo8-A_Inverse_resultat_scilab",inv(A));





