%(* ::Package:: *)

% Self-consistent DFT program for the helium atom.
% Solve using equidistant finite differences.
% Combine Hartree + Schroedinger. Self-consistent loop.

clear all;

format long
% Construct mesh with Nmax mesh points. 
% rmin and rmax are the end points of the mesh.
rmin = 0.00001; % first mesh point
rmax = 7; % last mesh point
Nmax = 200; % number of mesh points
h = (rmax - rmin)/(Nmax-1);

r=linspace(rmin,rmax,Nmax);

% hydrogen or helium?
q = 2; % helium
% q = 1; % hydrogen

% Convergence criteria
convcrit=1e-07;


% Construct first guess of charge density * r
for i = 2:Nmax-1,
   r2DensityInit(i) = % ## ## # fill in here ## ## #; % a hydrogen 1s orbital, suitably normalized.
end
r2DensityInit(1) = 0; r2DensityInit(Nmax) = 0;

% Check the normalization: q/4pi = integral (r^2*phi^2) = integral (r2Density)
Normr2DensityInit = trapz(r,r2DensityInit);

r2DensityOld = r2DensityInit;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% hartree potential %%%%%%%%%%%
%%%%%%%%%%%% CE1 %%%%%%%%%%%%%%%%%%%%%%%

% Form of problem: Ay = b. 

% Construct trigonal square matrix AH. AH has dimension (Nmax-2).

AH=-2*eye(Nmax-2,Nmax-2);
for i = 1:Nmax-3,
   AH(i,i+1) = 1;
end
for i = 2:Nmax-2,
   AH(i,i-1) = 1;
end

% Boundary conditions
BC0 = 0; BCNmax = q;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% start of self-consistent loop here: take the calculated density and use as input.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ediff=100;
niter=1;
while (ediff >=convcrit)   % iteratate until convergence
   
   clear bH;
   clear yH;
   
% Calculate b on the mesh, and shift the index 
for i = 2:Nmax-1,
     bH(i-1) = -h^2*4*pi*r2DensityOld (i)/r(i);
end

% Add the boundary conditions 
bH(1) = bH(1) - BC0;
bH(Nmax-2) = bH(Nmax-2) - BCNmax;
% Make b into a column vector
bH = bH';

% Solve Ay=b to get y.
yH = AH\bH;

% final result - attach the first and last points.
% rVH is the hartree potential times r.
rVH(1) = BC0; rVH(Nmax) = BCNmax;
for i = 1:Nmax-2,
   rVH(i+1) = yH(i);
end

rVHsave(niter,:) = rVH; % save to check convergence of the Hartree potential

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Exchange term %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:Nmax,
   DensityOld(i) = r2DensityOld (i)/(r(i)*r(i));
   ex(i) = % ## ## # fill in here ## ## #
   Vx(i) = % ## ## # fill in here ## ## #
%   Vx(i) = 0; ex(i) = 0; % test the effect of removing exchange
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Correlation term %%%%%%%%%%%%%%%%
%%%%%%% Thijssen %%%%%%%%%%%%%%%%%%%%%%
% Unpolarised
gamma = -0.1423; beta1 = 1.0529; beta2 = 0.3334; 
A = 0.0311; B = -0.048; C = 0.0020; D = -0.0116;
%%%%%%%%%%%%%
% Polarised
% gamma = -0.0843; beta1 = 1.3981; beta2 = 0.2611; 
% A = 0.01555; B = -0.0269; C = 0.0014; D = -0.0108;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:Nmax,
   DensityOld(i) = r2DensityOld (i)/(r(i)*r(i));
   if DensityOld(i) < 0.0001,
      Vc(i) = 0; ec(i) = 0;
   else  
      rs = % ## ## # fill in here ## ## #
      if rs >= 1,
         ec(i) = gamma/(1 + beta1*sqrt(rs) + beta2*rs); % rs > 1
         Vc(i) = ec(i)*(1+(7/6)*beta1*sqrt(rs)+beta2*rs)/(1+beta1*sqrt(rs)+beta2*rs); % rs > 1
      else
         ec(i) = A*log(rs) + B + C*rs*log(rs) + D*rs; % rs < 1
         Vc(i) = A*log(rs) + B - A/3 + (2/3)*C*rs*log(rs) + (2*D-C)*rs/3; % rs < 1
      end
   end
%   Vc(i) = 0; ec(i) = 0; % test the effect of removing correlation
   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Schroedinger equation %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Schroedinger equation with potential terms (Hartree, external, and exchange-correlation)
% Form of problem: Differential equation -1/2*y'' - (-VH(x) + q/x ) * y = lambda*y. lambda = eigenvalue. 
% q = charge (q=2 for helium, q=1 for hydrogen)
% VH = hartree potential
% more general form: -0.5y'' + Vy = lambda*y, V = potential.

% Construct tridiagonal matrix A. This square matrix has dimension Nmax-2.
for i = 1:Nmax-2,
   for j = 1:Nmax-2,
      AS(i,j) = 0;  
      AS(i,i) = -2 - 2*h^2*( rVH (i+1)/r(i+1) + Vx(i+1) + Vc(i+1) - q/r(i+1) ); % form -2-2h^2*V, where V is the potential 
%      AS(i,i) = -2 - 2*h^2*( - 1/r(i+1) ); % Simplified form of the Hamiltonian, for problem 2, hydrogen atom
   end
end
for i = 1:Nmax-3,
   AS(i,i+1) = 1;
end
for i = 2:Nmax-2,
   AS(i,i-1) = 1;
end

% rescale A so that we get an ordinary eigenvalue problem
AS = AS./(-2*h^2);

% Boundary conditions: The wave function phi is zero at first and last point, gives y(1) = y(Nmax) = 0,
% enabling the simple matrix eigenvalue formulation of the equation.
% Solve (A + eE)y = 0. E = unity matrix.   (See CE2)
lambda = eig(AS);
[lambdaMin, ilambdaMin] = min(lambda); % find the ground state
[v,d] = eig(AS); % column ilambdaMin of v should be the ground-state eigenvector.
ySMin = v(:,ilambdaMin);

% re-index
for i = 1:Nmax-2,
   rphi(i+1) = ySMin(i);
end

% set the end point values to the boundary condition values.
rphi(1) = 0; rphi(Nmax) = 0;

% Normalize the calculated (rphi)^2 solution to q/(4*pi)
r2Density = rphi.*rphi.*(q/(4*pi))./( trapz(r,rphi.*rphi) ); % This is the radial part of the density normalization.

% Copy over the new density. Mix.
mix = 0.5; % mixing parameter
r2DensityOld = % ## ## # fill in here ## ## #

% Calculate the total energy (Eq .32)
Etotal = q*lambdaMin - 8*pi*trapz(r, (1/q)*r2Density.*(0.5*rVH./r - ex - ec + Vx + Vc) )
Etot(niter,2) = Etotal; Etot(niter,1) = niter;
Eigenvalue(niter,2) = lambdaMin; Eigenvalue(niter,1) = niter;

  if(niter>=2) ediff=abs((Etot(niter,2)-Etot(niter-1,2))); end
  niter=niter+1;
end 
%%%%%%%%%%%%%%%%%%
% End of self-consistent loop here
%%%%%%%%%%%%%%%%%%

% analytical solutions for comparison: hydrogen r^2*R_ {10}^2, i.e., density*r^2
% for i = 1:Nmax, r2R2exact(i) = 4*pi*(r (i)^2)/pi*exp(-2*r(i)); end



figure(3)
clf
hold on
plot(r,r2Density,'k-')
xlabel('r (au)')
ylabel('r^2R^2 (au)')
box on

figure(4)
clf
plot(Etot(:,1),Etot(:,2),'k*-')
xlabel('Iteration number')
ylabel('Etot (Ha)')
box on

figure(5)
clf
plot(Eigenvalue(:,1),Eigenvalue(:,2),'r*-')
xlabel('Iteration number')
ylabel('\lambda _ {min}')
box on




