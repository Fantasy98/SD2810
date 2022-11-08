
% lwmain.m
%
% main program for aeroelastic analysis
%
% (c) 2004-2016 Dan Borglund <dodde@kth.se> and David Eller <dlr@kth.se>
% Note for Nov8 
clear all;

% setup geometry and structural properties
% number of finite elements requested should be a multiple of 3
nelem = 5;
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
rhop = 2000; % kg/m^3

% A guess of E and G
E = 22E9; % Gpa (from other group)
possion = 0.2;
G = E/2*(1+possion);
%G = E/2;

% set up linear constraints for clamped wing root
% Number of Degree of freedom
ndof = 3*nnodes;
B = eye(3,ndof);

npmass = 0 ;
dpm = zeros(npmass,3);
dpm(1,:) = 0;
dpm(2,:) =0;
% C : Null space of B  Chapter 6.8
%C = null(B);

% retrieve system matrices

[M1,K1,Z,Qip,f,CRv,CRd,s] = labwing(B, l, b, t, ba, mhinge, rhop, E, G, nelem, dpm);
fprintf("Offset s = %.2f m \n",s);
[V1,D1] = eig(K1,M1);
f1 = sqrt(abs(D1))./(2*pi);

% Retreive The result Again for trim the parameter
[M2,K2,Z,Qip,f,CRv,CRd,s] = labwing(B, 0.9*l, b, t, ba, mhinge, rhop, 0.8*E, G, nelem, dpm);
fprintf("Offset s = %.2f m \n",s);
[V2,D2] = eig(K2,M2);
f2 = sqrt(abs(D2))./(2*pi);

fprintf("The frequency of 1 is \n"); 
f1(1:5,1:5)
fprintf("The frequency of 2 is \n"); 
f2(1:5,1:5)