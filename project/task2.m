
% lwmain.m
% Task1 is for implementing the code for validation w
% Consider a semi-chord wing clamped at the root
% main program for aeroelastic analysis
%
% (c) 2004-2016 Dan Borglund <dodde@kth.se> and David Eller <dlr@kth.se>

clear all;

% setup geometry and structural properties
% number of finite elements requested should be a multiple of 3
nelem = 24;
nnodes = nelem + 1;

% lab wing dimensions and properties
# Semi-span
l =12.1/2; % m
% S = 7.41/2; % m^2
c = 0.64;
# According to the technic report, the main wing is assumed as ellipse with chord = 500mm
b = c/2 ; % m
# The alieron span is defined as c*E, E = 0.225 in W12C-JARA-0001.pdf P
ba = c*0.225/2; % m


# mass of hinge need to be figured out 2022/12/12
# Maybe use the same value as labwing is ok
mhinge = 28e-3;


% thickness is assumed to be 
t = 6e-4;
% t = c * .17;
# Unsure yet
% rhop = 1963.7; 
% m_wing = 80; %kg in W12C-JARA-0001.pdf P2
% rhop = 0.5*m_wing/(S*t);
# Density of carborn fiber -- load.pdf P6 , the results are quite close to previous 
rhop = 1760;



#Bending and torsional stiffness#
# Assume that the spar has the domiant contribution to both bending and torision
# The Shear Moduls has been given in load.pdf P9
G = 8600E6;
# According to its possion ratio v = 0.21 ~ 0.28, we estimate the Young's Moudles:
% E = 2*G*(1+0.28);
# According to spar.pdf P17
E = 23.9E9;

% definition matrix for discrete point masses to attach
% All attachment of contorl surface should be assigned as dpm
% The format is dpm(i,:) = [mass x_coord y_coord]
% 7 connection in total
% From  root to tip = 1:7
npmass = 7 ;
m1 = (40.33+6.39+2)/1000; %kg
m2 = (20.06+2*2)/1000; %kg
m3 = (40.33+2*6.39+2*2)/1000; %kg
x_coord = b -ba;
% x_coord = b-ba;
dpm = zeros(npmass,3);
dpm(1,:) = [m1 x_coord 0];
dpm(2,:) =[m2 x_coord 1*l/6];
dpm(3,:) =[m3 x_coord  2*l/6];
dpm(4,:) =[m2 x_coord  3*l/6];
dpm(5,:) =[m3 x_coord  4*l/6];
dpm(6,:) =[m2 x_coord 5*l/6];
dpm(7,:) =[m1 x_coord l];
% ....

% set up linear constraints for clamped wing root
% Number of Degree of freedom
ndof = 3*nnodes;
% Constrain of first 3 dof
B = [];

% retrieve system matrices

[M,K,Z,Qip,f,CRv,CRd] = nwing(B, l, b, t, ba, mhinge, rhop, E, G, nelem, dpm);



e1 = zeros(ndof,1);
e1(1:3:end) = 1;
e3 = zeros(ndof,1);
e3(3:3:end) = 1;


Q0 = Qip.Qtab(:,:,1);
# Test at ulitimate load corresponds to the diving speed 
u = 350/3.6;
g = 9.81;
# Load factor = 9 
nz = 9;
q = 0.5*Qip.rho* u * u;
# Consider the clamped wing condition as mentioned in structure test
% yn = 0.5*l*linspace(-1,1,nnodes);
% B = zeros(3,ndof);
% B(1,:) = e1'

% B(2,2:3:end) = 1;

% B(2,1:3:end) = yn;

% idof = 3*(nelem/2)+3;
% B(3,idof) = 1;
B(1,:) = e1;
B(2,2) = 1;
B(3,3) = 1;
Z = null(B);


LHS = [ Z' *(K-q*Q0)*Z  -q*Z'*Q0*e3
        q*e1'*Q0*Z       q*e1'*Q0*e3];

% Right hand side

RHS = [-nz*g*Z'*M*e1 
        g*e1'*M*e1];

x = LHS\RHS;
alfa0 =x(end) 
vhat = x(1:end-1);
vtot =  Z*vhat + alfa0*e3;


