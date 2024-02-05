function [D] = fdmatrix(x,diff_ord,accu_ord)

% Function to output a finite difference matrix D based on the function fdcoefs
%   x is the grid of nodes [x0 x1 x2...xN]
%   diff_ord is the differentiation order
%   accu_ord is the accuracy order

n = accu_ord+diff_ord-1;    %n is the number of points in the stencil
m = diff_ord;               %m is the differentiation order

npoints = length(x);        %Number of grid points
nedge = floor((n+1)/2);     %Number of end points requiring consideration 

D = zeros(npoints);

if mod(n,2)==1
    for i = 1:nedge
        D(i,1:(n+1)) = fdcoefs(m,n,x(1:(n+1)),x(i));
    end
    for i = (nedge+1):(npoints-nedge)
        D(i,(i-nedge):(i+nedge-1)) = fdcoefs(m,n,x((i-nedge):(i+nedge-1)),x(i));
    end
    for i = (npoints-nedge+1):(npoints)
        D(i,(end-n):end) = fdcoefs(m,n,x((end-n):end),x(i));
    end
else
    for i = 1:nedge
        D(i,1:(n+1)) = fdcoefs(m,n,x(1:(n+1)),x(i));
    end
    for i = (nedge+1):(npoints-nedge)
        D(i,(i-nedge):(i+nedge)) = fdcoefs(m,n,x((i-nedge):(i+nedge)),x(i));
    end
    for i = (npoints-nedge+1):(npoints)
        D(i,(end-n):end) = fdcoefs(m,n,x((end-n):end),x(i));
    end
end

%D = sparse(D);

end % end functioon fdmatrix


function coefs = fdcoefs(m,n,x,xi)

% Finite difference weights (Fornberg)
%   m: Differentiation order
%   n: size of the stencil (3,5,7,9)
%   x: coordinates of the points in the stencil
%   xi: coordinate of the evaluation point
% Number of points in the formula n+1: formal order n-m+1 (irregular points)

c1= 1;
c4= x(1)-xi;

c= zeros(n+1,m+1);
c(1,1)= 1;

for i=1:n-1;
    mn= min([i,m]);
    c2= 1;
    c5= c4;
    c4= x(i+1)-xi;
    for j=0:i-1;
        c3= x(i+1)-x(j+1);
        c2= c2*c3;
        for k= mn:-1:1;
            c(i+1,k+1)= c1*(k*c(i,k)-c5*c(i,k+1))/c2;
        end;
        c(i+1,1)= -c1*c5*c(i,1)/c2;
        for k=mn:-1:1;
            c(j+1,k+1)= (c4*c(j+1,k+1)-k*c(j+1,k))/c3;
        end;
        c(j+1,1)= c4*c(j+1,1)/c3;
    end;
    c1= c2;
end;

coefs= c(:,m+1)';

end % end function fdcoefs
