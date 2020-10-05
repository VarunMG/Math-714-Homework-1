%x_j goes from j = 0,......,N+2
%y_j goes from j = 0,......,N
%the matrix for -Delta_h will be of size (N+1)^2 x (N+1)^2
%%%%%%fine grid
format long;
N = 257;
dx = 1/(N+2);
dy = 1/N;
% e = ones(N+1,1);
% A = spdiags([e,-2*e,e], -1:1, N+1,N+1);
% B = spdiags([e,-2*e,e], -1:1, N+1,N+1);
% B(1,2) = 2; B(end,end-1) =2;
% A = (-1/(dx*dx))*A;
% B = (-1/(dy*dy))*B;
% Id = speye(N+1);
% Delh = kron(Id,A) + kron(B,Id);
% F = zeros((N+1)^2,1);
% count = 0;
% for i = 1:N+1:(N+1)^2
%     F(i) = f(count*dy)/(dx*dx);
%     count = count + 1;
% end
%u = Delh\F;
%u = jacobi_iter(Delh,F,20000);
u = jacobi_iter(N,50000);
[X,Y] = meshgrid(0:dx:1,0:dy:1);
%%%%% coarser grid
N2 = 65;
Newdx = 1/(N2+2);
Newdy = 1/N2;
[X2,Y2] = meshgrid(0:Newdx:1,0:Newdy:1);
u2 = jacobi_iter(N2,50000);
interpVals = interp2(X2,Y2,u2,X,Y,'linear'); %linear interpolant onto finer mesh for comparison
AbsErr = abs(interpVals-u);
disp('Max-norm Err:')
max(max(AbsErr))