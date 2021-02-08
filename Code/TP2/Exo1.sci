function [x]=solsup(U,b)
    n = size(U)(1); //Nombre lignes de matrice
    for i=n:-1:1
        if abs(U(i,i)) < %eps then
            error("Matrice U est non inversible");
        end
        sum_temp = 0;
        for j=i+1:n
            sum_temp = sum_temp + U(i,j)*x(j);
        end
        x(i) = 1/U(i,i)*(b(i)-sum_temp);
    end
endfunction
