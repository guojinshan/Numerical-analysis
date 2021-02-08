/*
Créateur: Jinshan GUO et Anais Debureaux
*/
exec(fullpath(pwd() + '\TP6.sci'),-1);

disp("===========Exo1===========");
t = [0 1 3 4 5.5 6];
tau_1 = t;
A1 = construct(t,tau_1);
disp("Resulat_Exo1(c)_A1=", A1);

tau_2 = [0.5 2 2.5 3.5 4.5 5.75];
A2 = construct(t,tau_2);
disp("Resulat_Exo1(d)_A2=", A2);

tau_3 = 0.3 * ((1:21) - 1);
A3 = construct(t,tau_3);
disp("Resulat_Exo1(e)_A3=", A3);

disp("===========Exo2===========");
    // Calculer la solution  equations normales
y = [1 1.5 1.25 0 0 1.5];
z = mcnorm(A2, y);
disp("Resulat_Exo2(a)_z=",z);
g_z = A2 * z;
    // Vérifier g_z
disp("Resulat_Exo2(b_i)_g_z=",g_z);
disp("Le résultat semble correcte car A est une matrice du rang plein");
    //  Tracer
plot(t,z, 'r');
scatter(tau_2, y, marker=4);
xgrid;
title("Application2",'fontsize',4);
legend('Parabole obtenue','Point coordonnées',1);
    // Calculer la solution  equations normales
y=[0 0.6 1.4 1.7 2.1 1.9 1.6 1.4 1.4 1 0.5 0.4 -0.2 -0.8 -0.5 0 0.4 1 1.6 1.7 1.2];
z = mcnorm(A3, y);
g_z = A3 * z;
//  Tracer
scf();
plot(t, z, 'r');
scatter(tau_3, y, marker=4);
xgrid;
title("Application3",'fontsize',4);
legend('Parabole obtenue','Point coordonnées',1);

condidtionnement =  cond(A3);
disp("Resulat_Exo2(x_ii)_cond=",condidtionnement);

disp("===========Exo3===========");
tol = 10^-3;
N = 100;
x0 = [1 2]'
[x, k] = newton(foncjac, tol, N, x0)
disp("Resulat_Exo3(x)=",x);
disp("Resulat_Exo3(k)=",k);

courbe(10)

disp("===========Exo4===========");
disp("On doit résoudre Rz=c1 avec QTy=(c1; d1)");
y=[0 0.6 1.4 1.7 2.1 1.9 1.6 1.4 1.4 1 0.5 0.4 0.2 0.8 0.5 0 0.4 1 1.6 1.7 1.2]
z = mcQR(A3,y);
disp("Resulat_Exo2_z=", z);
[Q,R]=qr(A3);
n = size(A3,2);
Rt=R(1:n,:);
disp("Resulat_Exo4(Rt_ii)_cond=",cond(Rt));
disp("Resulat_Exo4(M_ii)_cond=",cond(M));


