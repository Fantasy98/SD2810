
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
l =1.6; % m
b = 0.175; % m
ba = 0.03; % m
% ba = 0;
% measured from lab


% mhinge = ( 2*(40.33+6.39+2)+...
%             2*(20.06+2*2)+...
%             2*(40.33+2*6.39+2*2)+...
%             3*28.9)  % g
% mhinge = mhinge/1000; % kg 

% Here the mhinge is just mass of one rod 
% In the function labwing, it will be tripped
mhinge = 28e-3;

% mhinge = 0;
t = 0.004;%m
rhop = 1963.7; % Measured Density

% Measured E and G by viberation test
E = 31.5E9;
G = 5.52E9;

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
omega = sqrt(LAMBDA);


##fprintf("Eigen Frequencies are %.2f rad/s \n",LAMBDA);
% show modeshapes

##plotmode(V);

% stop here until the rest is implemented


% compute divergence speed
[udiv,zdiv] = divergence(K, Qip);
fprintf(1,'Divergence speed: %.2f m/s \n', udiv);
%compute reversal speed
[urev,zrev] = reversal(K, Qip, f, CRv, CRd);
fprintf(1,'Reversal speed: %.2f m/s \n', urev);

nmode = 3;
neig = 3
[Km,Mm,Zm,mQip]=ReduceDim(M,K,Qip,nmode);
[ucrit,pcrit,zcrit,pconv,uvec] = flutter(Mm,Km,mQip,neig,50);

fprintf(1,'Fluteer speed: %.2f m/s \n', ucrit);

return;

##% compute flutter speed
##[ucrit, pcrit, zcrit] = flutter(M,K,Qip);
##fcrit =
##fprintf(1,'Flutter speed: %.2f m/s \n',ucrit);
##fprintf(1,'Frequency of the critical mode: %.2f Hz \n',fcrit);
##
##% look at flutter mode shape
##vismode( ??? );
