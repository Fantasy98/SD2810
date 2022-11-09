
% lwmain.m
%
% main program for aeroelastic analysis
%
% (c) 2004-2016 Dan Borglund <dodde@kth.se> and David Eller <dlr@kth.se>

clear all;

% setup geometry and structural properties
% number of finite elements requested should be a multiple of 3
nelem = 12;

nnodes = nelem + 1;

% lab wing dimensions and properties
l =1.2; % m
b = 0.2; % m
ba = 0.03; % m
ba = 0;
% measured from lab
mhinge = (40.3 + 38.78 + 20.06 + 6.39 + 28.92 + 20.06 + 10.39 + 4.08 + 5.57)/1000; % kg
mhinge = 0;
t = 0.003;%m
rhop = 2000; % m^3/kg
rhop = 1900;
% A guess of E and G 
E = 25E9;
% Test Result
E = 20E9; % Gpa
% Assumed Possion Ratio
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
%% Number of Degree of freedom
ndof = 3*nnodes;
B = eye(3,ndof);
% Here we are testing condition without clamping
%B = [];
% retrieve system matrices
[M,K,Z,Qip,f,CRv,CRd,s] = labwing(B, l, b, t, ba, mhinge, rhop, E, G, nelem, dpm);

fprintf("Offset s = %.2f m \n",s);

% Compute the Inneria
I= (2*b*t^3)/12

% print free vibration frequencies
[V,LAMBDA] = eig(K,M);
Vhat = Z * V;
omega = sqrt(Z*LAMBDA)./(2*pi);
% fprintf("Eigen Frequencies are %.2f rad/s \n",omega(1:5,1:5));

% Compute beam mass per length
my = rhop * b*t;
% Constant Bn values for  3 lowest modes for free
muL = [4.730 7.853 (2*3+1)*pi/2];
f_bending = (muL.^2) * sqrt(E*I/(my*l^4)) ;

fprintf("Frequency of bending is %.2f rad/s \n", f_bending );
% show modeshapes
%% Plot the bending mode
vismode(Vhat(:,1))

% stop here until the rest is implemented
return;

