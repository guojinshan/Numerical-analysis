/*
Créateur: Jinshan GUO et Anais Debureaux
*/
exec(fullpath(pwd() + '\TP4.sci'),-1);

disp("===========Exo1-Factorisation de Cholesky===========");
A1 = [15 10 18 12;10 15 7 13;18 7 27 7;12 13 7 22];
[C]=cholesky(A1);
    //Exo1_a
disp("Resulat_Exo1(a)_C=",C);
disp("Si A1 == A1'':", A1 == A1');
disp("Car A1 est symétrique et tous les termes diagonaux de C sont postives, A1 est SDP");
    //Exo1_b
b1 = [53;72;26;97];
[x]=resolchol(A1,b1);
disp("Resulat_Exo1(b)_x=",x);

disp("===========Exo2-Matrice Hibert===========");
    //Exo2_a
H4 = hilbermat(4);
H10 = hilbermat(10);
    //Exo2_b
disp("Resulat_Exo2(b)_H4=",cond_norm_2(H4));
disp("Resulat_Exo2(b)_H10=",cond_norm_2(H10));
    //Exo2_c
disp("Resulat_Exo2(c)_H4_scilab=", cond(H4));
disp("Resulat_Exo2(c)_H10_scilab=", cond(H10));
disp("Donc on obtient le même résultat que celui de scilab");
    //Exo2_d
b = [1 1 1 1]';
[x]=resolchol(H4,b);
disp("Resulat_Exo2(d)_x=",x);
b_tild = [0.999 1.004 0.995 1.005]';
[x_tild]=resolchol(H4,b_tild);
disp("Resulat_Exo2(d)_x_tild=",x_tild);
disp("Resulat_Exo2(d)_EcartRelatifB=",norm(abs(b_tild-b),2)/norm(b,2));
disp("Resulat_Exo2(d)_EcartRelatifX=",norm(abs(x_tild-x),2)/norm(x,2));
    //Exo2_e
b1 = ones(10,1);
b1_tild = rand(b1,"uniform") * 10^-3;
[x1]=resolchol(H10,b1);
disp("Resulat_Exo2(e)_x=",x1);
[x1_tild]=resolchol(H10,b1_tild);
disp("Resulat_Exo2(e)_x_tild=",x1_tild);
disp("Resulat_Exo2(e)_EcartRelatifB=",norm(abs(b1_tild-b1),2)/norm(b1,2));
disp("Resulat_Exo2(e)_EcartRelatifX=",norm(abs(x1_tild-x1),2)/norm(x1,2));
disp("Resulat_Exo2(e)_x_scilab=",H10\b1);
disp("Resulat_Exo2(e)_x_tild_scilab=",H10\b1_tild);


disp("===========Exo3-Point fixe===========");
    //Question b
x = linspace(0,%pi,100);
[C]=pointfixe(cos,100,10^-12,1);
plot(x, cos(x),'b+');
plot(x, x,'r-');
plot(C,C,'go');



/*
disp("===========Exo4-Matrice creuse===========");
n1 = 10;
n2 = 50;
    //Cas A2
A2_n10 = Constuire_A2(n1);
disp("Resulat_Exo4_A2(n=10)=",A2_n10);
disp("Resulat_Exo4_A2_C(n=10)=",cholesky(A2_n10));
A2_n50 = Constuire_A2(n2);
disp("Resulat_Exo4_A2(n=50)=",A2_n50);
disp("Resulat_Exo4_A2_C(n=50)=",cholesky(A2_n50));   
    //Cas A3
A3_n10 = Constuire_A3(n1);
disp("Resulat_Exo4_A3(n=10)=",A3_n10);
disp("Resulat_Exo4_A3_C(n=10)=",cholesky(A3_n10));
A3_n50 = Constuire_A3(n2);
disp("Resulat_Exo4_A3(n=50)=",A3_n50);
disp("Resulat_Exo4_A3_C(n=50)=",cholesky(A3_n50));
*/


