/*
Code à vérifier le résultat dans un fichier .sce
*/ 


funcprot(0);
exec(fullpath(pwd() + '\Exo_pre.sci'),-1);
T = [0,1,3,4];
cc = [1 2 1;4 1 -1;1 -4 1];
clf();
[t1,z1] = trace(10,T,cc);
subplot(2,2,1);
plot(t1,z1);
xtitle("Figure1 N=10");

[t2,z2] = trace(100,T,cc);
subplot(2,2,2);
plot(t2,z2);
xtitle("Figure2 N=100");

[t3,z3] = trace(200,T,cc);
subplot(2,2,3);
plot(t3,z3);
xtitle("Figure3 N=200");

[t4,z4] = trace(400,T,cc);
subplot(2,2,4);
plot(t4,z4);
xtitle("Figure4 N=400");

