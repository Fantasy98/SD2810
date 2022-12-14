
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
ba = 0;
% ba = 0;
% ba = 0.2;
% measured from lab % Not correct need to be calculate again
mhinge = (40.3 + 38.78 + 20.06 + 6.39 + 28.92 + 20.06 + 10.39 + 4.08 + 5.57)/1000; % kg
mhinge = 0;
t = 0.003;%m

rhop = 1963.7; % Measured Density

% Measured E and G by viberation test
E = 25E9;
% Assumed Possion Ratio
E = E * 2.5;
possion = 0.21;
G = E/2*(1+possion);


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
[Modes,Pressure] = eig(A,K);
qDiv = 1/max(abs(diag(Pressure)));
fprintf("The divergence pressure is %.2f pa \n",qDiv);
air_rho = Qip.rho;
UDiv = sqrt(2*qDiv/air_rho);


fprintf("The divergence speed is %.2f m/s \n",UDiv);
