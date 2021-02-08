/*
Créateur: Jinshan GUO et Anais Debureaux
*/

function [z] = eulerexpl(f, y0, t0, T, N)
    p = length(y0); // y0 est un vecteur de Mp,1
    h = (T-t0)/N; // Calculer le pas
    z = zeros(p,N+1); 
    z(:,1) = y0; // Initialisation de z
    t = linspace(t0,T,N+1); // Constuire t
    for n=1:N
        z(:,n+1) = z(:,n) + h * f(t(n),z(:,n));
    end
endfunction

function [z] = eulerimpl(f, y0, t0, T, N)
    p = length(y0); // y0 est un vecteur de Mp,1
    h = (T-t0)/N; // Calculer le pas
    z = zeros(p,N+1); 
    z(:,1) = y0; // Initialisation de z
    t = linspace(t0,T,N+1); // Constuire t
    for n=1:N
        // résolution de  z(:,n+1) = z(:,n) + h * f(t(n+1),z(:,n+1));
        // z(:,n+1) = u est solution de l'équation : u = z(:,n) +  h * f(t(n+1),u)
        // g(u) = 0 ==> g(u) = -u + z(:,n) +  h * f(t(n+1),u)
        deff('[z] = g(u)','z = -u + z(:,n) +  h * f(t(n+1),u)');
        z(:,n+1) = fsolve(z(:,n) ,g); // resol g(u) = 0
    end
endfunction

// Initialisation (f, y0, t0, T, N)
deff('[dy] = f(t,y)','dy = y + t^2'); //y'(t) = y(t) + t^2
y0 = 0; //Condition initiale
t0 = 0;
T = 5;

// Verification de la solution exacte
//y'(t) = 2e^t - 2t - 2 = y(t) + t^2

function plot_figure(f, y0, t0, T, N)
    t = linspace(t0,T,N+1);
    // Solution euleur explicite
    z_eulerexpl =  eulerexpl(f, y0, t0, T, N);
    // Solution euleur implicite
    z_eulerimpl =  eulerimpl(f, y0, t0, T, N);
    // Solution exacte
    z_exacte = 2 * exp(t) - t^2 - 2*t - 2;
    // Solution donnée par ODE
    z_ode = ode(y0,t0,t,f);
    // Figure
    scf();
    plot2d(t, [z_eulerexpl', z_eulerimpl',z_exacte', z_ode'],leg='Euler Explicite@Euler Implicite@Exacte@ODE');
    
    // Figure Erreur
    err_eulerexpl = abs(z_eulerexpl - z_exacte);
    err_eulerimpl = abs(z_eulerimpl - z_exacte);
    err_ode = abs(z_ode - z_exacte);
    figure();
    plot2d(t,[err_eulerexpl',err_eulerimpl', err_ode'], leg='Erreur Euler Explicite@Erreur Euler Implicite@Erreur ODE');
endfunction

// Cas N=10
N1 = 10;
plot_figure(f, y0, t0, T, N1);

// Cas N=20
N2 = 20;
//plot_figure(f, y0, t0, T, N2);

// Cas N=50
N3 = 50;
//plot_figure(f, y0, t0, T, N3);

