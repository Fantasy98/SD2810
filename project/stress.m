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
t = 6e-3;
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
dpm = zeros(npmass,3);
% set up linear constraints for clamped wing root
% Number of Degree of freedom
ndof = 3*nnodes;

B = eye(3,ndof);
Z = null(B);

[M,K,Z,Qip,f,CRv,CRd] = nwing(B, l, b, t, ba, mhinge, rhop, E, G, nelem, dpm);

Ixx = 1.4117e+05 * 10^(-12);
I = Ixx;
Pload = -1;
P = zeros(ndof,1);

P(1:3:end-2) = Pload;
% P(end) = Pload;

P_hat = P' * Z;
v = K \ P_hat';


v = Z * v;
% plotmode(v(:,1));
delta_estimate = v(end-2);
fprintf("The FEM solution is %.5f \n",delta_estimate);
% Compute the Inneria

% Deformation
delta_theory = P(end-2)*l^3 /(3*E*I);

fprintf("The Analytical Solution is %.5f\n",delta_theory);


############## Plot bending moment and normal stress
ndof = length(v);

nnode = fix(ndof/3);

nshape = [1, nnode];
w = zeros(nshape);
w1 = zeros(nshape);
theta = zeros(nshape);

for k = 1:nnode
    w(k) = v(1+3*(k-1)); % nodal deflection
    w1(k) = v(2+3*(k-1)); % nodal deflection curvature    
    theta(k) = 180/pi * v(3+3*(k-1)); % nodal twist
end

we = zeros(nnode*2,1);
we(1:2:end) =  v(1:3:end); 
we(2:2:end) =  v(2:3:end-1);

yp = linspace(0.0, l, nnode);
ye = linspace(0,l,nelem);
le = yp(2)-yp(1);
w2 = diff(w1)/le;
t = c * .17;
Mx = -E * Ixx .* w2';
sigma = ye'.*Mx /Ixx;
plot(ye,Mx,"ro-","linewidth",1.5);

figure(2)
plot(ye,sigma,"bo-","linewidth",1.5);

