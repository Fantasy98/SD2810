
% Given the deflection to plot the distribution
% (c) 2004-2016 Dan Borglund <dodde@kth.se> and David Eller <dlr@kth.se>

clear all;

% setup geometry and structural properties
% number of finite elements requested should be a multiple of 3
nelem = 24;
nnodes = nelem + 1;

% lab wing dimensions and properties
l =1.6; % m
b = 0.175; % m
ba = 0.03 ;
% measured from lab % Not correct need to be calculate again
mhinge = 0;
t = 0.004;%m

rhop = 1950; % Measured Density

% Measured E and G by viberation test
E = 31.5E9;
% Assumed Possion Ratio
G = 5.52E9;


% definition matrix for discrete point masses to attach
npmass = 0 ;
dpm = zeros(npmass,3);
dpm(1,:) = 0;
dpm(2,:) =0;
% ....

% set up linear constraints for clamped wing root
% Number of Degree of freedom
ndof = 3*nnodes;
B = eye(3,ndof);

% retrieve system matrices

[M,K,Z,Qip,f,CRv,CRd,s] = labwing(B, l, b, t, ba, mhinge, rhop, E, G, nelem, dpm);
fprintf("Offset s = %.2f m \n",s);

[V,D] = eig(K,M);
Vhat = Z * V;
v = Vhat(:,1);

ndof = length(v);
nnode = fix(ndof/3);
nshape = [1, nnode];

% Extract bending and torsion dofs section by section
w = zeros(nshape);
t = zeros(nshape);
for k = 1:nnode
  w(k) = v(1+3*(k-1)); % nodal deflection
  t(k) = v(3+3*(k-1)); % nodal twist
end

