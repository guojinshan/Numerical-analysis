/*
TP7B Monome06: Jiawen ZHU
Ven. 8 jan
*/
//---------------Exo ptMilieu---------------
function [Y] = ptMilieu(a, t0, T, N, f)
    p = length(a);
    Y = zeros(p,N);
    Y(:,1) = a;
    h = T/(N-1);
    t = t0;
    for i=1:N-1
        Y(:,i+1) = Y(:,i) + h*f(t+h/2,Y(:,i)+h*f(t,Y(:,i))/2)
        t = t+h;
    end
endfunction

//---------------Ex fvdp--------------------
function out = fvdp(t,y)
    c = 0.4
    A = [0, 1;
        -1, c*(1-y(1)^2) ]
//    out = A * y
    out = [y(2);
            c*(1-y(1)^2)*y(2) - y(1)]
endfunction

//---------------euler-------------------
function [Y] = euler(a, t0, T, N, f)
    p = length(a);
    Y = zeros(p,N);
    Y(:,1) = a
    h = T/(N-1);
    t = t0;
    for i=1:N-1
        Y(:,i+1) = Y(:,i) + h*f(t,Y(:,i));
        t = t+h;
    end
endfunction
//---------------RK4---------------
function [X] = RK4(a, t0, T, N, f)
    p = length(a);
    X = zeros(p,N);
    X(:,1) = a
    h = T/(N-1);
    t = linspace(t0,t0+T,N-1)
    K = zeros(p,4);
    for i = 1:N-1
         K(:,1) = f(t(i),X(:,i))
         K(:,2) = f(t(i)+h/2,X(:,i)+h*K(:,1)/2)
         K(:,3) = f(t(i)+h/2,X(:,i)+h*K(:,2)/2)
         K(:,4) = f(t(i)+h,X(:,i)+h*K(:,3))
         
         X(:,i+1) = X(:,i) + h*(K(:,1) + 2* (K(:,2)+ K(:,3)) + K(:,4))/6
    end
endfunction

//---------------VDP---------------
function tracevdp(a,t0,T,Neuler,Nptmil,NRK4,Node)
    teuler = linspace(t0,t0+T,Neuler) //abscisses
    tptmil = linspace(t0,t0+T,Nptmil) //abscisses
    tRK4 = linspace(t0,t0+T,NRK4) //abscisses
    tode = linspace(t0,t0+T,Node) //abscisses

    Yeuler = euler(a, t0, T, Neuler, fvdp)
    Yptmilieu = ptMilieu(a, t0, T, Nptmil, fvdp)
    YRK4 = RK4(a, t0, T, NRK4, fvdp)
    Yode = ode(a,t0,tode,fvdp)

    plot(Yeuler(1,:),Yeuler(2,:),'-.g',Yptmilieu(1,:),Yptmilieu(2,:),'--k',YRK4(1,:),YRK4(2,:),':b',Yode(1,:),Yode(2,:),'-.r');
    legend("Euler","PtMil","RK4","ODE");
    axe=get("current_axes");
//    axe.data_bounds=[t0,min(Y(1,$));T,max(Y(1,$))];
endfunction
//---------------TV---------------
function tab=TV(a,t0,T,f)
    tab = zeros(4,4);
    over = ['   N','     Euler','     PTM','       RK4']
    for i=1:4
        N = 10^i
        t = linspace(t0,t0+T,N) //abscisses
        Yeuler = euler(a, t0, T, N, f)
        Yptmilieu = ptMilieu(a, t0, T, N, f)
        YRK4 = RK4(a, t0, T, N, f)
        
        tab(i,1) = N
        tab(i,2) = Yeuler(1,$)
        tab(i,3) = Yptmilieu(1,$)
        tab(i,4) = YRK4(1,$)
    end
    disp(over)
    disp(tab)
endfunction

//---------------TE---------------
function [tab]=TE(a,t0,T,f)
    tab = zeros(4,4);
    over = ['   N','     Euler','     PTM','       RK4']
    for i=1:4
        N = 10^i
        t = linspace(t0,t0+T,N) //abscisses
        Yeuler = euler(a, t0, T, N, f)
        Yptmilieu = ptMilieu(a, t0, T, N, f)
        YRK4 = RK4(a, t0, T, N, f)
        
        val = solution(t($))
        tab(i,1) = N
        tab(i,2) = abs(val -Yeuler(1,$))
        tab(i,3) = abs(val -Yptmilieu(1,$))
        tab(i,4) = abs(val -YRK4(1,$))
    end
    disp(over)
    disp(tab)
endfunction
function out = solution(t)
    out = 3/2*exp(-t) -1/2*cos(t)+1/2*sin(t)
endfunction

function out = f(t,y)
    out = -y+sin(t)
endfunction
//---------------plot---------------
function ordre = traceErreur(tabErreur)
    logVal = zeros(length(log(tabErreur(:,1))),4);
    logVal(:,1) = log(tabErreur(:,1));

    logVal(:,2) = -log(tabErreur(:,2));
    logVal(:,3) = -log(tabErreur(:,3));
    logVal(:,4) = -log(tabErreur(:,4));
    figure();
    plot2d(logVal(:,1),[logVal(:,2),logVal(:,3),logVal(:,4)], leg='Erreur Euler Explicite@Erreur Point Milieu@Erreur RK4');
    plot(logVal(:,1),logVal(:,2),'ro');
    plot(logVal(:,1),logVal(:,3),'y+');
    plot(logVal(:,1),logVal(:,4), 'b>');
    h1 = legend(['Erreur Euler Explicite','Erreur Point Milieu','Erreur RK4','MC Euler Explicite','MC Point Milieu','MC RK4'],-1);
    ordre = zeros(3,1)
     xlabel("log10(N)");
    ylabel("−log10(Erreur)");
    for i = 2:4 
       [a,b] = reglin(logVal(:,1)',logVal(:,i)')
        ordre(i-1)=a;
    end
endfunction


