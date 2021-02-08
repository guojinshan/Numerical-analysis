function [U,e]=trisup(A,b)
    n = size(A)(1); //Nombre de ligne de matrice
    for k=1:n-1
        if abs(A(k,k)) < %eps then
            error('Matrice est non inversible pendant élimination.');
        else
            for i = (k+1):n
                coef = A(i,k)/A(k,k);
                b(i) = b(i) - coef*b(k);
                A(i,k) = 0;
                for j=(k+1):n
                    A(i,j) = A(i,j) - coef*A(k,j);
                end
            end
        end      
    end
    if abs(A(n,n)) < %eps then
        error("Le dernier pivot est null, Matrice non inversible pendant l''élimination.");
    end
    U = A;
    e = b;
endfunction

