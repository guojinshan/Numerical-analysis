 function [L,U] = LU(A)
     n = size(A)(1) //Nombre de ligne de matrice
     //Initialisation de L et U
     L = eye(n,n);
     U = zeros(n,n);
     for k=1:n
         for j=k:n
             temp = 0;
             for p=1:k-1
                temp =  temp + L(k,p) * U(p,j);
             end
             U(k,j) = A(k,j) - temp;
         end
         if abs(U(k,k)) < %eps then
            error('Pivot nul');
         end
         for i=k+1:n
             temp1 = 0;
             for p=1:k-1
                temp1 = temp1 + L(i,p) * U(p,k);
             end
              L(i,k) = 1/U(k,k) * (A(i,k) -temp1);
         end
     end
 endfunction
