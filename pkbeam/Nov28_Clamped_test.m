
% Nov28 
% Concept of Mass matrix
% Primilary Experimet for validating the clamped wing in wind tunnel
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
npmass = 7
m1 = (40.33+6.39+2)/1000; %kg
m2 = (20.06+2*2)/1000; %kg
m3 = (40.33+2*6.39+2*2)/1000; %kg
x_coord = (280-175)/1000;
dpm = zeros(npmass,3);
dpm(1,:) = [m1 x_coord 0];
dpm(2,:) =[m2 x_coord 27/100];
dpm(3,:) =[m3 x_coord  53/100];
dpm(4,:) =[m2 x_coord  80/100];
dpm(5,:) =[m2 x_coord  106/100];
dpm(6,:) =[m2 x_coord 133/100];
dpm(7,:) =[m2 x_coord 160/100];
% ....

% set up linear constraints for clamped wing root
% Number of Degree of freedom
ndof = 3*nnodes;
% Constrain of first 3 dof
B = eye(3,ndof);

% retrieve system matrices

[M,K,Z,Qip,f,CRv,CRd,s] = labwing(B, l, b, t, ba, mhinge, rhop, E, G, nelem, dpm);
fprintf("Offset s = %.2f m \n",s);

[V,LAMBDA] = eig(K,M);
Vhat = Z * V;
omega = diag(sqrt(LAMBDA)./(2*pi));
freq = real(omega(1:20));

for i = 0:4
    fprintf(" \n Frequency-%.2d:  \n %.2frad/s \n %.2frad/s \n %.2frad/s \n",i+1,freq(1+3*i),...
                                                                                    freq(2+3*i),...
                                                                                    freq(3+3*i))

    % fprintf(" 2 + 3N Bending frequencies is %.2f rad/s \n",freq(2+3*i))

    % fprintf(" 3 + 3N Torsional frequencies is %.2f rad/s \n",freq(3+3*i))
end


% % Computing the Wing mass 
% % we constrained mass, the mass of wing is calculated as: 
% e1 = zeros(ndof-3,1);
% e1(1:3:end) = 1; 
% mass_wing = e1' * M * e1;
% fprintf("mass of wing %.2f m \n",mass_wing);

% % Make sure the mass is constrained. when get the mass matrix each time

% % For Flutter: Plot u VS Imp as well !!

% % The result turns to be less than what calculated by labwing

% % Now let's use free constrain 
% B = [];
% [M_free,K,Z,Qip,f,CRv,CRd,s] = labwing(B, l, b, t, ba, mhinge, rhop, E, G, nelem, dpm);
% e1 = zeros(ndof,1);
% e1(1:3:end) = 1; 
% mass_wing_free = e1' * M_free * e1;
% fprintf("mass of wing free constrain %.2f m \n",mass_wing_free);

% return;


% % compute divergence speed
% [udiv,zdiv] = divergence(K, Qip);
% fprintf(1,'Divergence speed: %.2f m/s \n', udiv);
% %compute reversal speed
% [urev,zrev] = reversal(K, Qip, f, CRv, CRd);
% fprintf(1,'Reversal speed: %.2f m/s \n', urev);


% ##% compute flutter speed
% [ucrit, pcrit, zcrit] = flutter(M,K,Qip);
% fcrit = pcrit;
% fprintf(1,'Flutter speed: %.2f m/s \n',ucrit);
% fprintf(1,'Frequency of the critical mode: %.2f Hz \n',fcrit);
% ##
% ##% look at flutter mode shape
% vismode( zcrit );
