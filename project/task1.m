
% lwmain.m
% Task1 is for implementing the code for validation w
% Consider a semi-chord wing clamped at the root
% main program for aeroelastic analysis
%
% (c) 2004-2016 Dan Borglund <dodde@kth.se> and David Eller <dlr@kth.se>

clear all;

% setup geometry and structural properties
% number of finite elements requested should be a multiple of 3
nelem = 9;
nnodes = nelem + 1;

% lab wing dimensions and properties
# Semi-span
l =12.1/2; % m
S = 7.41/2; % m^2
c = 0.64;
# According to the technic report, the main wing is assumed as ellipse with chord = 500mm
b = c/2 ; % m
# The alieron span is defined as c*E, E = 0.225 in W12C-JARA-0001.pdf P
ba = c*0.225/2; % m


# mass of hinge need to be figured out 2022/12/12
# Maybe use the same value as labwing is ok
mhinge = 6*28E-3;


% thickness is assumed to be 
t = 0.947 * 10^(-3);
t = 6E-3;

% t = 1e-3
% t = 1E-3;
# Unsure yet
% rhop = 1963.7; 
% m_wing = 80; %kg in W12C-JARA-0001.pdf P2
% rhop = 0.5*m_wing/(S*t);
# Density of carborn fiber -- load.pdf P6 , the results are quite close to previous 
rhop = 1760;
% rhop = 1200


#Bending and torsional stiffness#
# Assume that the spar has the domiant contribution to both bending and torision
# The Shear Moduls has been given in load.pdf P9
G = 8600E6;
# According to its possion ratio v = 0.21 ~ 0.28, we estimate the Young's Moudles:
% E = 2*G*(1+0.28);
# According to spar.pdf P17
E = 23.9E9;



npmass = 7 ;
m1 = 6*(40.33+6.39+2)/1000; %kg
m2 = 6*(20.06+2*2)/1000; %kg
m3 = 6*(40.33+2*6.39+2*2)/1000; %kg
x_coord = (b-ba);
dpm = zeros(npmass,3);
dpm(1,:) = [m1 x_coord 0];
dpm(2,:) =[m2 x_coord l/6];
dpm(3,:) =[m3 x_coord  2*l/6];
dpm(4,:) =[m2 x_coord  3*l/6];
dpm(5,:) =[m3 x_coord  4*l/6];
dpm(6,:) =[m2 x_coord 5*l/6];
dpm(7,:) =[m1 x_coord l];


% set up linear constraints for clamped wing root
% dpm = zeros(npmass,3);
% Number of Degree of freedom
ndof = 3*nnodes;
% Constrain of first 3 dof
% retrieve system matrices
B = [];


[M,K,Z,Qip,f,CRv,CRd] = nwing(B, l, b, t, ba, mhinge, rhop, E, G, nelem, dpm);


e1 = zeros(ndof,1);
e1(1:3:end) = 1;
e3 = zeros(ndof,1);
e3(3:3:end) = 1;

mtot = e1' * M * e1;


B = eye(3,ndof);
% B(1,:) = e1;
Z= null(B);



# Simulate distribution force
P = zeros(ndof,1);
P_dis = 10727/l;
% P(1:3:end-2) = 0.8*P_dis;
P(1:3:nnodes*0.5*3-2) = 0.8*P_dis;
P(nnodes*0.5*3+1:3:end-2) = 0.6 * P_dis;
P(3:3:end) = -1095;
P_hat = P' * Z;
v = (Z' * K * Z) \ P_hat';

% The mode to be plotted by plotmode() 
v_mode  = Z * v;
% plotmode(v_mode(:,1));
delta_estimate = v(end-2);
fprintf("The FEM solution is %.5f \n",delta_estimate);
% Compute the Inneria

% Deformation
delta_theory = P(end-2)*l^3 /(3*E*I);

fprintf("The Analytical Solution is %.5f\n",delta_theory);

plot_stress(Z*v,l,b,t,E,G,c)
print -djpg -r0 task1s.jpg
return;
