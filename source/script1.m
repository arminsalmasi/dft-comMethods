% Calculate the Hartree potential from the hydrogen atom charge density.
% Form of problem: Ay = b, where y contains y(2) to y(npts-1), and b
% contains b(2) to b(npts-1).


clear all;

% Construct mesh with npts mesh points. rmin and rmax are the end points of the mesh.


rmin = 0;
rmax = 5;
npts = 100;

r=linspace(rmin,rmax,npts);
h = %##### insert correct expression for the step size h #####

% Construct trigonal matrix A. This square matrix has dimension npts-2.
% Zero the matrix

A=-2*eye(npts-2,npts-2);          %Diagonal
for i = 1:npts-3,                 %Upper diagonal
   A(i,i+1) = 1;
end
for i = 2:npts-2,                 %Lower diagonal
   A(i,i-1) = 1;
end

% Boundary conditions on b
q = 1;
bc0 = %##### insert correct boundary condition #### ; 
bcnpts = %##### insert correct boundary condition ####;

% Calculate b on the mesh, and shift the index 
for i = 2:npts-1,
   b(i-1) = %#### insert correct expression for b ##### % one 1s orbital in the density

end

% Add the boundary conditions 
b(1) = b(1) - bc0;                          %y(0)
b(npts-2) = b(npts-2) - bcnpts;             %y(rmax)

% Make b into a column vector
b = b';

% Solve Ay=b
y = A\b;

% final result
yfinal(1) = bc0; yfinal(npts) = bcnpts;
for i = 1:npts-2,
   yfinal(i+1) = y(i);
end

% get the potential from yfinal
rVH = yfinal;

% control. Is this expression really the correct one??
for i = 1:npts,
   rVHtheor(i) = 1-exp(-2*r(i))*(r(i)+1);
end

figure(1)
clf
hold on
plot(r,rVH,'k-')
plot(r,rVHtheor,'r-')
legend('Calc','Analytical',1)
xlabel('r(au)')
ylabel('rV_H (au)')
box on

