function [x]=solinf(L,b)
    n = size(L)(1); //Nombre lignes de matrice
    for i=1:n
        if abs(L(i,i)) < %eps then
            error("Matrice L est non inversible");
        end
        sum_temp = 0;
        for j=1:i-1
            sum_temp = sum_temp + L(i,j)*x(j);
        end
        x(i) = 1/L(i,i)*(b(i)-sum_temp);
    end
endfunction
