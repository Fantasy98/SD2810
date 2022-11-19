
% Given the deflection to plot the distribution
% (c) 2004-2016 Dan Borglund <dodde@kth.se> and David Eller <dlr@kth.se>

clear all;

% setup geometry and structural properties
% number of finite elements requested should be a multiple of 3
nelem = 20;
nnodes = nelem + 1;

% lab wing dimensions and properties
l =1.6; % m
b = 0.175; % m
ba = 0;
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
%Let's start with a given load and compare to analytical result!
Pload = -1;
P = zeros(ndof,1);
% Add a point load 

P(end-2) = Pload;
% Add a torque
P(end) = Pload;

P_hat = Z' * P;
v = K \ P_hat;

ndof = length(v);
nnode = fix(ndof/3);
nshape = [1, nnode];

% Extract bending and torsion dofs section by section
w = zeros(nshape);
theta = zeros(nshape);
for k = 1:nnode
  w(k) = v(1+3*(k-1)); % nodal deflection
  w2(k) = v(2+3*(k-1)); % nodal deflection curvature
  
  theta(k) = 180/pi * v(3+3*(k-1)); % nodal twist
end

% Compute normal stress along span
% Moment = -E * I * w'',
% w'' is the curvature of beam
% Then the normal stress can be derived as sigma = M * z/I
% z = 0.5 * t if we only consider the outer part
yp = linspace(0.0, l, nnode);
Ixx= (2*b*t^3)/12;


% delta_theory = P(end-2)*l^3 /(3*E*Ixx)

Mx = -E * Ixx .* w2 .*yp;
sigma = Mx .*0.5*t./Ixx;

% Compute shear stress along span
% theta = T * L / (J*G)
beta = 1/3;
alpha = 1/3;
Tx = theta*beta * 2*b * t^3 * G ./yp;
tau = Tx./(alpha * (2*b) * t^2);
clf;
figure(1);
subplot(2,1,1);
hold on 
plot(yp,sigma(end:-1:1),"ro-","linewidth",1.5);
title("Normal Stress Distribution","fontsize",15)
ylabel("Normal Stress (pa)","fontsize",8)


subplot(2,1,2);
plot(yp,Mx(end:-1:1),"bs-","linewidth",1.5)
title("Bending Moment Distribution","fontsize",15)
ylabel("Moment (N*m)","fontsize",8)


figure(2);
subplot(2,1,1);
hold on 
plot(yp,tau(end:-1:1),"ro-","linewidth",1.5);
title("Shear Stress Distribution","fontsize",15)
ylabel("Normal Stress (pa)","fontsize",8)


subplot(2,1,2);
plot(yp,Tx(end:-1:1),"bs-","linewidth",1.5)
title("Shear Force Distribution","fontsize",15)
ylabel("Shear Force (N*m)","fontsize",8)

figure(3)
plotmode(v);

% Minimum Safety Factor = Ultimate_stress/ max_Allowed_stress
% In our case, since the material property is similar to Aluminum, we can assume 
% Assume the ultimate stress is identical to alumn's 
Ultimate = 90e6;
factor_normal = Ultimate / max(abs(sigma));
factor_shear = Ultimate / max(abs(tau(2:end)));
fprintf("Minimum Safety factor of normal stress is %.2f\n", factor_normal)
fprintf("Minimum Safety factor of shear stress is %.2f\n", factor_shear)