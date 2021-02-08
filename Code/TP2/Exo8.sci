function [B]=inverse(A)
    n = size(A)(1); //Nombre lignes de matrice
    for k = 1:n
        pivot_max = A(k,k);
        nombre_ligne = k;
        for i=k:n
            if A(i,k) > pivot_max then
                nombre_ligne = i;
                pivot_max = A(i,k);
            end
        end
        if Est_Inversible == %F then
            error("Matrice non inversible");
        else
            //Echanger la ligne i et la ligne k
            temp = A(k,:);
            A(k,:) = A(nombre_ligne,:);
            A(nombre_ligne,:) = temp;
            
        end
    end
endfunction
