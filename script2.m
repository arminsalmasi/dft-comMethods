% Solve Schrodinger equation (for hydrogen atom) using finite differences.
% Differential equation: -1/2*y'' - 1/x * y = ey. e = eigenvalue. 
% Solve using equidistant finite differences.

%% #Armin: I followed the instruction here. The output plot shows analytical 
%          and calcualted values are very close
%%

clear variables; clc; close all;

% Construct mesh with npts mesh points. rmin and rmax are the end points of the mesh.
rmin = 0.00000001;
rmax = 20;
npts = 400;

r = linspace(rmin,rmax,npts);
h = (rmin - rmax)/npts ;% h is the step size
%% #Armin: It may be more efficient to use toeplitz function to make a toeplitz
% Construct tridiagonal matrix A. This square matrix has dimension npts-2.
% Zero the matrix
A=-2*eye(npts-2,npts-2);          %Diagonal
B = zeros(1,npts-2);
B(1)= -2; B(2) = 1;
A = toeplitz(B);

% rescale A so that we get an ordinary eigenvalue problem
A = A./(-2*h^2);

% Solve (A + eE)y = 0. E = unity matrix. 
e1 = eig(A);
[emin,iemin] = min(e1); % find the ground state
[v,d] = eig(A); % column iemin of v should be the ground-state eigenvector.
%% #Armin: v(:,iemin) is the right eigen vectore here which is the solution! 
ymin = v(:,iemin); 

% re-index
for i = 1:npts-2,
   yminfinal(i+1) = ymin(i);
end
%% #Armin: after applying boundry conditions, yminfinal = r*phi
% Boundary conditions. Wave function is zero at first and last point, y(1) = y(npts) = 0,
% enabling the simple matrix eigenvalue formulation of the equation.
% set the end point values to the boundary condition values.
yminfinal(1) = 0; yminfinal(npts) = 0;

% and then we want the charge density:  rho = |phi|^2, yminfinal = r * phi
phi = yminfinal./r;
rho = abs(phi) .* abs(phi);
r2R2Calc = yminfinal .* yminfinal;

% check the norm
rhoNorm = trapz(r,rho);
yNorm = trapz(r,yminfinal);

% analytical solutions for comparison
% hydrogen r^2*R_{10}^2, i.e., density*r^2

for i = 1:npts,
   r2R2Anal(i) = (r(i)^2)/pi*exp(-2*r(i));
end
r2R2Norm = trapz(r,r2R2Anal); % = 1 as expected when the factor 4pi is included. 


% Normalize the calculated r2R2 solution to 1/4*pi
r2R2Calc = r2R2Calc./(trapz(r,r2R2Calc)*4*pi);



figure(2)
clf
hold on
plot(r,r2R2Calc,'k-')
plot(r,r2R2Anal,'r-')
plot(r,rho,'b-')

legend('Calc','Analytical','\rho=n(r)')
xlabel('r(au)')
ylabel('r^2R^2 (au) & \rho = n(r) = \phi^2')
box on
