exec('C:\Users\Administrator\Desktop\MT09\Code\TP1\Exercice1.sci',-1);
N = 60;
x1 = suite(N,0,1);
x2 = suite(N,1/3,1/12);
clf;
subplot(1,2,1);plot(1:60,x1','ro');xtitle("Itérés de la suite quand a=0 et b=1","X","Y");
subplot(1,2,2);plot(1:60,x2','g+');xtitle("Itérés de la suite quand a=1/3 et b=1/12","X","Y");





