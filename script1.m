% Calculate the Hartree potential from the hydrogen atom charge density.
% Form of problem: Ay = b, where y contains y(2) to y(npts-1), and b
% contains b(2) to b(npts-1).
%% #Armin: I followed the instruction here. The output plot shows analytical 
%          and calcualted values are very close 
%%

clear variables; clc; close all;

% Construct mesh with npts mesh points. rmin and rmax are the end points of the mesh.
rmin = 0;
rmax = 5;
npts = 100;

r = linspace(rmin,rmax,npts);
%% #Armin: I am not sure about abs() term. It may be completely unncessary 
h = abs(rmin - rmax)/npts ;% insert correct expression for the step size h
%% #Armin: It may be more efficient to use toeplitz function to make a toeplitz
% Construct trigonal matrix A. This square matrix has dimension npts-2.
% Zero the matrix
A=-2*eye(npts-2,npts-2);          %Diagonal
B = zeros(1,npts-2);
B(1)= -2; B(2) = 1;
A = toeplitz(B);
%           for i = 1:npts-3,                 %Upper diagonal
%              A(i,i+1) = 1;
%           end
%           for i = 2:npts-2,                 %Lower diagonal
%              A(i,i-1) = 1;
%           end
%% #Armin: Boundary conditions on b
q = 1; %hydrogen % n umber of electrons
bc0 = 0; %##### insert correct boundary condition #### ; by solving charge from potential inside a large sphear
bcnpts = q; %##### insert correct boundary condition ####;
%% #Armin: Calculate b on the mesh, and shift the index % Poisson equation
%          U(r) = rVH(r) gives q = 4.Pi.int(r^2 .n(r). dr) [normalized electron density]
%          [d^2U(r)/dr^2 = -4.Pi.r.n(r)] = first term & [n(r)= exp(-2.r)/Pi] gives
%          [d^2U(r)/dr^2 = -4.Pi.r.(exp(-2.r)/Pi]  
%          [h^2 comes from finite diference method
%          y''(xi)=(yi+1?2yi+yi?1)/h^2]
%          [b is d^2U(r)/dr^2 = U''(r) = {rVH(r)}'']
%          aim is solving U''(r) to find U(r) = rVH = VH = the Hartree
%          potential so numerically b(i-1) = h^2 * (-4.pi.r.n(r))
for i = 2:npts-1
   b(i-1) = -4*r(i)*exp(-2*r(i))*h^2 ; %one 1s orbital in the density !n = exp(-2r)/pi given in the text
end
%% 
% Add the boundary conditions 
b(1) = b(1) - bc0;                          %y(0)
b(npts-2) = b(npts-2) - bcnpts;             %y(rmax)

% Make b into a column vector
b = b';

% Solve Ay=b 
y = A\b;


% final result
yfinal(1) = bc0; yfinal(npts) = bcnpts;
for i = 1:npts-2
   yfinal(i+1) = y(i);
end

% get the potential from yfinal
rVH = yfinal;

% control. Is this expression really the correct one??
%for i = 1:npts
   rVHtheor(:) = 1-exp(-2*r(:)).*(r(:)+1);
%end

VH = rVH ./r;

figure(1)
clf
hold on
plot(r,rVH,'k-')
plot(r,rVHtheor,'r-')
plot(r,VH,'b-')
legend('Calc','Analytical', 'V_H')
xlabel('r(au)')
ylabel('rV_H (au) =U(r) & V_H')
box on

