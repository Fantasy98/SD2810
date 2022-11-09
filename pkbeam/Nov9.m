
% Nov 9 Divergence Analysis
clear all;

% setup geometry and structural properties
% number of finite elements requested should be a multiple of 3
nelem = 10;
nnodes = nelem + 1;

% lab wing dimensions and properties
l =1.6; % m
b = 0.175; % m
% b = 0.2 ;
ba = 0.03; % m
% ba = 0;
% ba = 0.2;
% measured from lab % Not correct need to be calculate again
mhinge = (40.3 + 38.78 + 20.06 + 6.39 + 28.92 + 20.06 + 10.39 + 4.08 + 5.57)/1000; % kg
mhinge = 0;
t = 0.003;%m
rhop = 2000; % kg/m^3

% A guess of E and G
E = 22E9; % Gpa (from other group)
possion = 0.3;
G = E/2*(1+possion);
%G = E/2;

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

% C : Null space of B  Chapter 6.8
% C = null(B);

% retrieve system matrices

[M,K,Z,Qip,f,CRv,CRd,s] = labwing(B, l, b, t, ba, mhinge, rhop, E, G, nelem, dpm);
fprintf("Offset s = %.2f m \n",s);


% P the constrain added in the DOF
%
Pload = -1;
P = zeros(ndof,1);
P(end-2) = Pload;
P_hat = Z' * P;
v = K \ P_hat;
delta_estimate = v(end-2)

% Compute the Inneria
I= (2*b*t^3)/12

% Deformation
delta_theory = P(end-2)*l^3 /(3*E*I)

% print free vibration frequencies
[V,LAMBDA] = eig(K,M);
Vhat = Z * V;
omega = sqrt(abs(diag(LAMBDA)))./(2*pi);

% Check the struct Qip i.e dictionary in python
% fieldnames(Qip)
% Check the reduced frequency
% Qip.ktab
% For divergence == steady state use ktab = 0

%Qip.Qtab
% Aerodynamic force
% check the first , which corrsponding to steady state
% a = Qip.Qtab(:,:,1)
% Bending not generate aerodynamic force, 
% Twist will give the aerodynamic force

% spy() spartisiy of the data
% spy(a) . do not have data at n & 2n because bending will not bring aerodynamic force 
% spy(a) @ 3n has value because it due to the twist
% It is due to the strip theory which assume that lift is only related each section's lift and momentum
% For flat plate Cl_alpha = 2*pi

A = Qip.Qtab(:,:,1);
% Solve the divergence 
[Shape,Pressure] = eig(A,K);
air_rho = Qip.rho;
% Find the max eigenvalue and 
% the inverse is the divergence pressure
maxvalue = max(abs(diag(Pressure)));
qDiv = 1 ./ maxvalue
VDiv = sqrt(2*qDiv/air_rho)
