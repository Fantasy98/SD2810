
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
S = 7.41/2; % m^2
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
x_coord = (280-175)/1000;
% x_coord = b-ba;
dpm = zeros(npmass,3);
dpm(1,:) = [m1 x_coord 0];
dpm(2,:) =[m2 x_coord 27/100];
dpm(3,:) =[m3 x_coord  53/100];
dpm(4,:) =[m2 x_coord  80/100];
dpm(5,:) =[m3 x_coord  106/100];
dpm(6,:) =[m2 x_coord 133/100];
dpm(7,:) =[m1 x_coord 160/100];
% ....

% set up linear constraints for clamped wing root
% Number of Degree of freedom
ndof = 3*nnodes;
% Constrain of first 3 dof
B = [];
B=eye(3,ndof);  
% retrieve system matrices

[M,K,Z,Qip,f,CRv,CRd] = nwing(B, l, b, t, ba, mhinge, rhop, E, G, nelem, dpm);


% [M,K,Z,Qip,f,CRv,CRd] = labwing(B, l, b, t, ba, mhinge, rhop, E, G, nelem, dpm);






% compute divergence speed
[udiv,zdiv] = divergence(K, Qip);
fprintf(1,'Divergence speed: %.2f m/s \n', udiv);
fprintf(1,'1.15 times of maximum speed: %.2f m/s \n', 1.15*350/3.6);
fprintf("Divergence speed in reference report %.2f m/s \n", 670/3.6);
%compute reversal speed
% [urev,zrev] = reversal(K, Qip, f, CRv, CRd);
% fprintf(1,'Reversal speed: %.2f m/s \n', urev);



return;


nmode = 3;
neig = 3
[Km,Mm,Zm,mQip]=ReduceDim(M,K,Qip,nmode);
[ucrit,pcrit,zcrit,pconv,uvec] = flutter(Mm,Km,mQip,neig,50);

fprintf(1,'Fluteer speed: %.2f m/s \n', ucrit);
