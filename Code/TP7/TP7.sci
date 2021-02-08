/*
Créateur: Jinshan GUO et Anais Debureaux
*/
//===========Exo1 - Point Milieu ===========
function [z] = pointmilieu(a, t0, T, N, f)
    p = length(a); // a est un vecteur de Mp,1
    h = (T-t0)/N; // Calculer le pas
    z = zeros(p,N+1); // Initialisation de z
    zz = a; 
    z(:,1) = zz; // Stocker zz dans la première colonne de z
    t = linspace(t0,T,N+1); // Constuire t
    for n=1:N
        k1 = f(t(n),zz);
        k2 = f(t(n)+h/2, zz + h/2 * k1);
        zz = zz + h * k2;
        z(:,n+1) = zz;
    end
endfunction

//===========Exo2===========
//Transformation de l'équation différentielle d'ordre 2 vers celle d'ordre 1
//x(t) = [x1(t);x2(t)] = [u(t);u'(t)]
//x(t)' =  [u'(t);u''(t)] = [x2(t);c*(1-x1(t)^2)*x2(t)-x1(t)]
//==> x(t)' = F(t,x(t))   'équation différentielle d'ordre 1'
// avec condition initiale x(t0) = [u(t0);u'(t0)] = [α; β]

function [Y] = vdp(t, X) // X=x(t) ==> vdp ne depend pas de t
    X1 =  X(1);
    X2 =  X(2);
    c = 0.4; // On fixe c = 0.4
    Y = [X2;c * (1-X1^2) * X2 - X1];
endfunction

function [z] = eulerexpl(a, t0, T, N, f)
    p = length(a); // a est un vecteur de Mp,1
    h = (T-t0)/N; // Calculer le pas
    z = zeros(p,N+1); 
    z(:,1) = a; // Initialisation de z
    t = linspace(t0,T,N+1); // Constuire t
    for n=1:N
        z(:,n+1) = z(:,n) + h * f(t(n),z(:,n));
    end
endfunction

//===========Exo4 -  Runge Kutta d’ordre 4 ===========
function [z] = RK4(a, t0, T, N, f)
    p = length(a); // a est un vecteur de Mp,1
    h = (T-t0)/N; // Calculer le pas
    z = zeros(p,N+1); // Initialisation de z
    zz = a; 
    z(:,1) = zz; // Stocker zz dans la première colonne de z
    t = linspace(t0,T,N+1); // Constuire t
    for n=1:N
        k0 = f(t(n),zz);
        k1 = f(t(n)+h/2, zz + h/2 * k0);
        k2 = f(t(n)+h/2, zz + h/2 * k1);
        k3 = f(t(n)+h, zz + h * k2);
        zz = zz + h/6 * (k0 + 2*k1 + 2 * k2 + k3);
        z(:,n+1) = zz; // Stocker zz dans la n+1 ère colonne de z
    end
endfunction

function tracevdp(a, t0, T, Neul, Nptmil, Node, Nrk4)
   deff('[dx] = f(t,x)' ,'dx = vdp(t, x)'); // Définir f
   scf();
   // Point milieu
   if  Nptmil <> 0 then
       z_pointmilieu =  pointmilieu(a, t0, T, Nptmil, f);
       //disp(z_pointmilieu);
       plot2d(z_pointmilieu(1,:), z_pointmilieu(2,:),-8); // Diagame de phase
   end
  // Euler explicite
  if  Neul <> 0 then
       z_eulerexpl =  eulerexpl(a, t0, T, Neul, f)
       //disp(z_eulerexpl);
       plot2d(z_eulerexpl(1,:), z_eulerexpl(2,:),3); // Diagame de phase
   end
   // ODE
   if  Node <> 0 then
       theta = linspace(t0,T,Node+1); // Constuire theta
       z_ode = ode(a, t0, theta, f);
       //disp(z_ode);
       plot2d(z_ode(1,:), z_ode(2,:),-1); // Diagame de phase
   end
   // RK4
   if  Nrk4 <> 0 then
       z_RK4 =  RK4(a, t0, T, Nrk4, f);
       //disp(z_RK4);
       plot2d(z_RK4(1,:), z_RK4(2,:),12); // Diagame de phase
   end
   xlabel("x1 - u(t)");
   ylabel("x2 - u''(t)");
   title("Portrait de Phase d''équation de Van Der Pol");
   h1 = legend(['Point milieu';'Euler explicite';'ODE';'RK4'],-1);
endfunction

//===========Exo5===========
function [TV, TE] = compar(a, t0, T)
    // Initialisation
    m = length(a); // a est un vecteur de Mm,1
    TV = zeros(m,3);
    deff('[du] = f(t,u)' ,'du = -u + sin(t)'); //f(t, u) = -u + sin(t)
    N = [10;100;1000;10000];
    Erreur = zeros(length(N),3);
    Ordre = zeros(length(N),3);
    
    for i=1:length(N)
        // Solution euleur explicite & Point milieu & RK4
        z_eulerexpl =  eulerexpl(a, t0, T, N(i), f);
        z_pointmilieu =  pointmilieu(a, t0, T, N(i), f);
        z_RK4 =  RK4(a, t0, T, N(i), f);
        TV(:,1) = z_eulerexpl(:,$);
        TV(:,2) = z_pointmilieu(:,$);
        TV(:,3) = z_RK4(:,$);
        
        // Solution exacte
        t = linspace(t0,T,N(i)+1); // Constuire t
        z_exacte = 3/2 * exp(-t) - 1/2 * cos(t) + 1/2 * sin(t);
        
        // Erreur
        TE = zeros(m,3);
        TE(:,1) = abs(z_exacte($) - z_eulerexpl(1,$));
        TE(:,2) = abs(z_exacte($) - z_pointmilieu(1,$));
        TE(:,3) = abs(z_exacte($) - z_RK4(1,$));
        
        // Stocker les erreurs
        Erreur(i,:) = TE(1,:);
    end
    Erreur = -log10(Erreur);
    disp('Erreur',Erreur)
    // Figure
    figure();
    plot2d(log10(N),[Erreur(:,1),Erreur(:,2),Erreur(:,3)]);
    // Approche chacune de ces courbes en moindres carrés
    [a_eulerexpl, b_eulerexpl] = reglin(log10(N), Erreur(:,1));
    [a_pointmilieu, b_pointmilieu] = reglin(log10(N), Erreur(:,2));
    [a_RK4, b_RK4] = reglin(log10(N), Erreur(:,3));
    plot(log10(N),Erreur(:,1), 'ro');
    plot(log10(N),Erreur(:,2), 'y+');
    plot(log10(N),Erreur(:,3), 'b>');
    h1 = legend(['Erreur Euler Explicite','Erreur Point Milieu','Erreur RK4','MC Euler Explicite','MC Point Milieu','MC RK4'],-1);
    xlabel("log10(N)");
    ylabel("−log10(Erreur)");
    
    disp("Resulat_Exo5(b)_TV=", TV);
    disp("Resulat_Exo5(b)_TE=", TE);
    
    //Commenter les valeurs des pentes de ces droites
    disp("Resulat_Exo5(d)_Pente_Eulerexpl=", a_eulerexpl);
    disp("Resulat_Exo5(d)_Pente_Pointmilieu=", a_pointmilieu);
    disp("Resulat_Exo5(d)_Pente_RK4=", a_RK4);
endfunction
