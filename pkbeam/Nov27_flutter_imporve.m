% Nov 23 Improve the flutter speed by adding discrete mass on the wing
% mhinge has been updated
% pk_bounds and pk_bisect

% Updated derivation of flutter speed
clear all;

% setup geometry and structural properties
% number of finite elements requested should be a multiple of 3
nelem = 15;
nnodes = nelem + 1;
l =1.6; % m
b = 0.175; % m
ba = 0.03 ;

% Here the mhinge is just mass of one rod 
% In the function labwing, it will be tripped
mhinge = 28e-3;
t = 0.004;%m

rhop = 1950; % Measured Density

% Measured E and G by viberation test
E = 31.5E9;
% Assumed Possion Ratio
G = 5.52E9;


% definition matrix for discrete point masses to attach
npmass = 12 ;

m1 = (40.33+6.39+2)/1000; %kg
m2 = (20.06+2*2)/1000; %kg
m3 = (40.33+2*6.39+2*2)/1000; %kg
% According to the photo
x_coord = (280-175)/1000;

% extra mass we can use to improve flutter speed
m_extra = [50 100 150]/1000;

% Fixed mass
dpm = zeros(npmass,3);
dpm(1,:) = [m1 x_coord 0];
dpm(2,:) =[m2 x_coord 27/100];
dpm(3,:) =[m3 x_coord  53/100];
dpm(4,:) =[m2 x_coord  80/100];
dpm(5,:) =[m2 x_coord  106/100];
dpm(6,:) =[m2 x_coord 133/100];
dpm(7,:) =[m2 x_coord 160/100];

% To improve flutter speed
% Adding extra discrete mass on the wing 
% Note: There is no point for adding extra mass at tip and root! 

% There is a solution but may not be the optimal for improving the flutter speed
x_coord_extra = x_coord+0.06;

dpm(8,:) = [3*m_extra(3) -x_coord_extra 4*l/5];
dpm(11,:) = [3*m_extra(3) -x_coord_extra 5*l/6];

% Test before/after adding hingmass 
% In theory, since hinges are after the elastic axis, the flutter speed should be decreased
% dpm = zeros(npmass,3);

% Consider discrete mass of hinges, the flutter speed is 15.7
% Without mass of hinges, the flutter speed is 16.0


% set up linear constraints for clamped wing root
% Number of Degree of freedom
ndof = 3*nnodes;
B = eye(3,ndof);

% retrieve system matrices

[M,K,Z,Qip,f,CRv,CRd,s] = labwing(B, l, b, t, ba, mhinge, rhop, E, G, nelem, dpm);
fprintf("Offset s = %.2f m \n",s);


[ucrit,pcrit,zcrit] = flutter(M,K,Qip);

fprintf("\nFlutter Speed  is %.2f m/s\n",ucrit);
fprintf("\nFlutter Frequency  is %.2f rad/s\n",pcrit);

% Visualize the mode
vismode(zcrit);


[urev,zrev] = reversal(K,Qip,f,CRv,CRd);
fprintf("\nReversal Speed  is %.2f m/s\n",urev);

[udiv,zdiv] = divergence(K,Qip);
fprintf("\nDivergence Speed  is %.2f m/s\n",udiv);




