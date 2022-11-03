
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
l =1.2; % m
b = 0.2; % m
ba = 0.03; % m
% measured from lab
mhinge = (40.3 + 38.78 + 20.06 + 6.39 + 28.92 + 20.06 + 10.39 + 4.08 + 5.57)/1000; % kg
%mhinge = 0;
t = 0.003;%m
rhop = 200; % m^3/kg

% A guess of E and G
E = 20E9; % Gpa
possion = 0.19;
G = E/2*(1+possion)
G = E/2;

% definition matrix for discrete point masses to attach
npmass = 0 ;
dpm = zeros(npmass,3);
dpm(1,:) = 0;
dpm(2,:) =0;
% ....

% set up linear constraints for clamped wing root
%% Number of Degree of freedom
ndof = 3*nnodes;
B = [] ;
B = eye(3,ndof);

% C : Null space of B  Chapter 6.8
C = null(B);

% retrieve system matrices

[M,K,Z,Qip,f,CRv,CRd,s] = labwing(B, l, b, t, ba, mhinge, rhop, E, G, nelem, dpm);
fprintf("Offset s = %.2f m \n",s);

% P the constrain added in the DOF
%
P = [0 0 0 1 0 0]';
##K\ P(4:6)

% Compute the Inneria
I= (2*b*t^3)/12
% Deformation
delta = P(4)*l^3 /(3*E*I)

% print free vibration frequencies
[V,LAMBDA] = eig(K,M);
omega = sqrt(LAMBDA);
##fprintf("Eigen Frequencies are %.2f rad/s \n",LAMBDA);
% show modeshapes

##plotmode(V);

% stop here until the rest is implemented
return;

##% compute divergence speed
##[udiv,zdiv] = divergence(K, Qip);
##fprintf(1,'Divergence speed: %.2f m/s \n', udiv);
##
##% compute reversal speed
##[urev,zrev] = reversal(K, Qip, f, CRv, CRd);
##fprintf(1,'Reversal speed: %.2f m/s \n', urev);
##
##% compute flutter speed
##[ucrit, pcrit, zcrit] = flutter(M,K,Qip);
##fcrit =
##fprintf(1,'Flutter speed: %.2f m/s \n',ucrit);
##fprintf(1,'Frequency of the critical mode: %.2f Hz \n',fcrit);
##
##% look at flutter mode shape
##vismode( ??? );
