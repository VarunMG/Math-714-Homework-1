function u = jacobi_iter_test(N,tol)
% [n,~] = size(A);
% D = diag(diag(A));
% N = D-A;
% f = inv(D)*b;
% IterMat = inv(D)*N;
% if max(eigs(IterMat)) > 1
%     error('spectral radius too large for convergence')
% end
% x = zeros(n,1);
% for i=1:iters
%    x =  IterMat*x + f;
% end
dx = (1/(N+2));
dy = (1/N);
[X,Y] = meshgrid(0:dx:1,0:dy:1);
exact = @(x,y) (-(cosh(2*pi)./sinh(2*pi)).*sinh(2*pi.*x)+cosh(2*pi.*x)).*cos(2*pi.*y);
exactVals = exact(X,Y);
dx2 = (1/(N+2))*(1/(N+2));
dy2 = (1/N)*(1/N);
const = 1/((2/dx2) + (2/dy2));
f = @(y) cos(2*pi*y);
leftVals = f(linspace(0,1,N+1));
u  = zeros(N+1,N+2);
u = [leftVals.' u];
counter = 0;
exactVec = reshape(exactVals.',(N+3)*(N+1),1);
while norm(exactVec-reshape(u.',(N+3)*(N+1),1)) > tol
    counter = counter + 1;
    uNew  = zeros(N+1,N+2);
    uNew = [leftVals.' uNew];
    for i=1:N+1
        for j=2:N+2
            if i==1
                uNew(i,j) = const*(2*u(i+1,j)/dy2 + u(i,j+1)/dx2 + u(i,j-1)/dx2);
            elseif i==N+1
                uNew(i,j) = const*(2*u(i-1,j)/dy2 + u(i,j+1)/dx2 + u(i,j-1)/dx2);
            else
                uNew(i,j) = const*(u(i+1,j)/dy2 + u(i-1,j)/dy2 + u(i,j+1)/dx2 + u(i,j-1)/dx2);
            end
        end
    end
    u = uNew;
end
counter
end


