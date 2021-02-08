/*
TPB Monome06: Jiawen ZHU
Ven. 8 Jan.
*/
exec(fullpath(pwd() + '\tp7.sci'),-1);
a = [2; -2]
t0 = 0
T = 15
Nptmil = 100;
Neuler = 100;
Node = 100;
NRK4 = 100;
tracevdp(a,t0,T,Neuler,Nptmil,NRK4,Node);
TV(a,t0,T,f);
f = TE(a,t0,T,f);
disp(traceErreur(f));
