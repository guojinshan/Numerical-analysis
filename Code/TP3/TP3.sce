/*
Créateur: Jinshan GUO et Anais Debureaux
*/

exec(fullpath(pwd() + '\TP3.sci'),-1);
//===========Exo1-Richtmayer===========
a = [1 2 1 2];
b = [2 3 1 1 1];
c = [3 1 2 1];
u = [-1 0 0 3 1];
[x]=rich(a,b,c,u);
disp("Resulat_Exo1_x=",x);

//===========Exo2===========
T = [1 3 4 5 5 6]';
t = [-5 1.5 4.5 6 1 20];
disp("Resulat_Exo2_(a)=");
for k=1:size(t)(2);
    [i]=place(T,t(k));
    if i then
        mprintf('\tLors que t=%.1f, i=%i\n',t(k),i);
    else
    end
end

disp("Resulat_Exo2_(b)=");
T1 = [1 3 4.5 5 6]';
cc = [1 0 1 0;5 0 -8/9 0;3 0 16 0;7 0 -8 0];
t1 = [3 5];
for k=1:length(t1);
    [z] = calcg(t1(k), T1, cc);
    mprintf('\tLors que t=%.1f, z=%.1f\n',t1(k),z);
end

trace(10,T1,cc);
trace(200,T1,cc);


//===========Exo3-Spline Cubique===========
T2 = [1 3 4.5 5 6]';
y = [1 5 3 7 -1]'
[d]=cald(T2,y)
disp("Resulat_Exo3a_(d)=",d);

[cc1] = calcoef(T2, y)
disp("Resulate_Exo3b_(cc)=",cc1);

scf();
trace(10,T2,cc1);
trace(200,T2,cc1);









