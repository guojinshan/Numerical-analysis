/*
    Algorithme qui détermine la valeur de l'indice i, pour t donné, 
    à placer dans un vector ordonné T.
*/
function [i]=place(T,t)
    n = length(T);
    if t<T(1)| t>T(n) then
        error("La valeur de t est hors de portée");
    end
    if t == T(n) then
        i = n-1; 
        return;
    end
    i_min = 1;
    i_max = n;
    while i_max - i_min > 1
        i_mil = floor((i_max + i_min) / 2);
        if t >= T(i_mil) then
           i_min = i_mil;
        else
            i_max = i_mil;
        end
    end
    i = i_min;
endfunction

/*
    Tracer la courbe d'équation
*/
function [t,z]=trace(N,T,cc)
    n_T = length(T);
    //Vérifier le nombre de colonne de matrice est correct
    n_cc = size(cc)(2);
    if n_cc <> 3 then
        error("Le nombre de colonne de matrice est fausse");
    end
    //Construire le vecteur t
    t = linspace(T(1),T(n_T),N);
    z = zeros(1,N);
    for i=1:N
        [indice] = place(T,t(i));
        z(i) = cc(indice,1) + cc(indice,2)*t(i) + cc(indice,3)*t(i)^2;
    end    
endfunction







