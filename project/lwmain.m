
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
l =12.1/2; % m
S = 7.41; % m^2
c = S/l;
# According to the technic report, the main wing is assumed as ellipse with chord = 500mm
b = c/2 ; % m
ba = (c - 2*0.5); % m


# mass of hinge need to be figured out 2022/12/12
mhinge = 28e-3;

# thickness is assumed to be 0.6mm in loads.pdf-P12
% t = c*0.017%m
t = 0.6e-3
# Unsure yet
rhop = 1963.7; 
% rhop = 80/(S*t)

% Bending and torsional stiffness
% In loads.pdf-P12 GK_tot = 158.6*c^3 
% where c is in mm,
% Torsion is in Angle/mm in RAD
GK = 158.6 * (c*1000)^3 * 10^9
G = 8600E6
E = GK/G

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
dpm = zeros(npmass,3);
dpm(1,:) = [m1 x_coord 0];
dpm(2,:) =[m2 x_coord 27/100];
dpm(3,:) =[m3 x_coord  53/100];
dpm(4,:) =[m2 x_coord  80/100];
dpm(5,:) =[m3 x_coord  106/100];
dpm(6,:) =[m2 x_coord 133/100];
dpm(7,:) =[m1 x_coord 160/100];
% ....
% dpm = zeros(npmass,3);
% set up linear constraints for clamped wing root
% Number of Degree of freedom
ndof = 3*nnodes;
% Constrain of first 3 dof
B = eye(3,ndof);
B = [];

% retrieve system matrices

[M,K,Z,Qip,f,CRv,CRd,s] = labwing(B, l, b, t, ba, mhinge, rhop, E, G, nelem, dpm);
fprintf("Offset s = %.2f m \n",s);


% compute divergence speed
[udiv,zdiv] = divergence(K, Qip);
fprintf(1,'Divergence speed: %.2f m/s \n', udiv);
%compute reversal speed
% [urev,zrev] = reversal(K, Qip, f, CRv, CRd);
% fprintf(1,'Reversal speed: %.2f m/s \n', urev);
return;


nmode = 3;
neig = 3
[Km,Mm,Zm,mQip]=ReduceDim(M,K,Qip,nmode);
[ucrit,pcrit,zcrit,pconv,uvec] = flutter(Mm,Km,mQip,neig,50);

fprintf(1,'Fluteer speed: %.2f m/s \n', ucrit);

