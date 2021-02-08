/*
Créateur: Jinshan GUO et Anais Debureaux
*/
exec(fullpath(pwd() + '\TP7.sci'),-1);

a = [2;-2];
t0 = 0;
T = 15;
Nptmil = 100;
Neul = 100;
Node = 100;
Nrk4 = 100;
tracevdp(a, t0, T, Neul, Nptmil, Node, Nrk4);
//tracevdp(a, t0, T, 1000, 1000, 1000, 1000);

//===========Exo5 ===========
T = 10;
[TV, TE] = compar(a, t0, T);
