function u = jacobi_iter(N,iters)
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
dx2 = (1/(N+2))*(1/(N+2));
dy2 = (1/N)*(1/N);
const = 1/((2/dx2) + (2/dy2));
f = @(y) cos(2*pi*y);
leftVals = f(linspace(0,1,N+1));
u  = zeros(N+1,N+2);
u = [leftVals.' u];
for k=1:iters
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

