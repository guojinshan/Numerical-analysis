// Interpolation_example.sce
function y=f(x)
        y = sin(x);
    endfunction
    
    N=30
    xj = linspace(0, 2*%pi, N+1); // row vector of size (N+1)
    xj = xj'; // make the vector vertical
    // Assembling the van der Monde matrix A ...
    A = zeros(N+1,N+1); // define the size
    for j=1:N+1
        A(:,j) = xj.^(j-1); // vectorized operation on columns
    end;
    //  then computes the polynomial coefficients
    // by solving the van der Monde system
    F = f(xj); // the result is a column vector
    a = A \ F; // solving the linear system, get the coefs
// Defines a function that computed the
    // polynomial interpolation
    function y=polinterpolation(x,a)
    // interpolation at position x with coefs in vector a
    l = length(a);
    y = 0*x; // initialize y (at size of x)
    for i=1:l
        y = y + a(i)*x.^(i-1);
    end;
    endfunction
    // Then plot the function and the interpolation polynomial
    // in the interval [0,2*pi] (for example)
    xd = linspace(0,2*%pi, 200);
    plot(xd, f(xd), '-r'); // in red with a line
    plot(xd, polinterpolation(xd,a), '.-b'); // in blue, dotted line
    plot(xj, F, 'o'); // location of interpolation points;
    plot(xd, f(xd), '-r', 'LineWidth', 2); // in red with a line
    // done.
