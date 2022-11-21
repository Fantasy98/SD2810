
% lwmain.m
%
% main program for aeroelastic analysis
%
% (c) 2004-2016 Dan Borglund <dodde@kth.se> and David Eller <dlr@kth.se>

clear all;

% setup geometry and structural properties
% number of finite elements requested should be a multiple of 3
nelem = 1;
nnodes = nelem + 1;

% lab wing dimensions and properties
l =1.6; % m
b = 0.175; % m
ba = 0.03; % m
ba = 0;
% measured from lab
mhinge = 0;
t = 0.004;%m
% rhop = 1963.7; % Measured Density
rhop = 2000;
% Measured E and G by viberation test
E = 25E9;
% Assumed Possion Ratio
% E = E * 2.5;
% possion = 0.21;
% G = E/2*(1+possion);
G = E;


% definition matrix for discrete point masses to attach
npmass = 0 ;
dpm = zeros(npmass,3);
dpm(1,:) = 0;
dpm(2,:) =0;
% ....

% set up linear constraints for clamped wing root
%% Number of Degree of freedom
ndof = 3*nnodes;
B = eye(3,ndof);

% C : Null space of B  Chapter 6.8


% retrieve system matrices

[M,K,Z,Qip,f,CRv,CRd,s] = labwing(B, l, b, t, ba, mhinge, rhop, E, G, nelem, dpm);
fprintf("Offset s = %.2f m \n",s);


% P the constrain added in the DOF
%
Pload = -1;
P = zeros(ndof,1);
P(end-2) = Pload;
P_hat = P' * Z;
v = K \ P_hat';

% The mode to be plotted by plotmode() 
v_mode  = Z * v;
plotmode(v_mode(:,1));
delta_estimate = v(end-2);
fprintf("The FEM solution is %.5f \n",delta_estimate);
% Compute the Inneria
I= (2*b*t^3)/12;

% Deformation
delta_theory = P(end-2)*l^3 /(3*E*I);

fprintf("The Analytical Solution is %.5f\n",delta_theory);
% print free vibration frequencies
[V,LAMBDA] = eig(K,M);
Vhat = Z * V;
omega = sqrt(LAMBDA);


