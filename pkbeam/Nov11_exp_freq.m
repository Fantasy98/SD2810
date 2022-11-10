
% For test the frequency with analytical eigenfrequncy
% Free End Beam
clear all;

% setup geometry and structural properties
% number of finite elements requested should be a multiple of 3
nelem = 20;

nnodes = nelem + 1;

% lab wing dimensions and properties
l =1.6; % m
b = 0.175; % m
% Here we do not consider the aileron 
%and mass of hinge since the test will be apply on the flat plate
ba = 0;
mhinge = 0;
t = 0.003;%m
% The estimate data from David
% which should be adjust after experiment !
%################################
rhop = 2000; % m^3/kg
rhop = 1900;
% A guess of E and G 
E = 25E9;
% Test Result
E = 22E9; % Gpa
% Assumed Possion Ratio
possion = 0.3;
G = E/2*(1+possion);
%################################

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
% No constrain == Free 
B = [];
% Here we are testing condition without clamping

% retrieve system matrices
[M,K,Z,Qip,f,CRv,CRd,s] = labwing(B, l, b, t, ba, mhinge, rhop, E, G, nelem, dpm);

fprintf("Offset s = %.2f m \n",s);

% Compute the Inneria
I= (2*b*t^3)/12

% print free vibration frequencies
[V,LAMBDA] = eig(K,M);
Vhat = Z * V;
omega = diag(sqrt(LAMBDA)./(2*pi));
fprintf("First Three bending frequencies from eigenvalue is %.2f rad/s \n",omega(1:4:16))
fprintf("First Three bending frequencies from eigenvalue is %.2f rad/s \n",omega(2:4:16))
fprintf("First Three bending frequencies from eigenvalue is %.2f rad/s \n",omega(4:4:16))
% Theoritical Formula for frequncies

fprintf("Analytical Solution of Frequencies: \n")
% 1. Bending frequency: 
% Compute beam mass per length
my = rhop * 2*b*t;
% Constant Bn values for  3 lowest modes for free
uL = [0 4.730 7.853];
r = 3:10;
url = (2.*r+1/2).*pi;
UL = [uL url];
f_bending = ((UL).^2)*sqrt(E*I/(my*l^4));
fprintf("Frequency of bending is %.2f rad/s \n", f_bending );

return;

