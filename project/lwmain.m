
% lwmain.m
%
% main program for aeroelastic analysis
%
% (c) 2004-2016 Dan Borglund <dodde@kth.se> and David Eller <dlr@kth.se>

clear all;

% setup geometry and structural properties
% number of finite elements requested should be a multiple of 3
nelem = 24;
nnodes = nelem + 1;

% lab wing dimensions and properties
l =12.1; % m
b = 0.5*(7.41/12.1); % m
ba = 0.114/2; % m
% ba = 0;
% measured from lab
mhinge = ( 2*(40.33+6.39+2)+...
            2*(20.06+2*2)+...
            2*(40.33+2*6.39+2*2)+...
            3*28.9); % g
mhinge = mhinge/1000; % kg 

% mhinge = 0;
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
% Constrain of first 3 dof
B = eye(3,ndof);


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
omega = diag(sqrt(LAMBDA)/(2*pi));


##fprintf("Eigen Frequencies are %.2f rad/s \n",LAMBDA);
% show modeshapes

##plotmode(V);

% stop here until the rest is implemented


% compute divergence speed
[udiv,zdiv] = divergence(K, Qip);
fprintf(1,'Divergence speed: %.2f m/s \n', udiv);
return;
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
