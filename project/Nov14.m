
% Nov 14 : Aileron Reversal
clear all;

% setup geometry and structural properties
% number of finite elements requested should be a multiple of 3
nelem = 10;
nnodes = nelem + 1;

% lab wing dimensions and properties
l =1.6; % m
b = 0.175; % m
% b = 0.2 ;
ba = 0; % m
% measured from lab % Not correct need to be calculate again
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

% retrieve system matrices

[M,K,Z,Qip,f,CRv,CRd,s] = labwing(B, l, b, t, ba, mhinge, rhop, E, G, nelem, dpm);
fprintf("Offset s = %.2f m \n",s);


A = Qip.Qtab(:,:,1);
[Modes,Pressure] = eig(A,K);
qDiv = 1/max(abs(diag(Pressure)));
fprintf("The divergence pressure is %.2f pa \n",qDiv);
air_rho = Qip.rho;
UDiv = sqrt(2*qDiv/air_rho);


fprintf("The divergence speed is %.2f m/s \n",UDiv);


