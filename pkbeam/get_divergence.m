
% Programme for computing divergence speed
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
mhinge = ( 2*(40.33+6.39+2)+...
            2*(20.06+2*2)+...
            2*(40.33+2*6.39+2*2)+...
            3*28.9)  % g
mhinge = mhinge/1000; % kg 

t = 0.003;%m

rhop = 1963.7; % Measured Density
fprintf("Reference Mass is 494g \n");

E = 25E9;
% Assumed Possion Ratio
E = E * 2.5;
possion = 0.21;
G = E/2*(1+possion);

% definition matrix for discrete point masses to attach
npmass = 0 ;
dpm = zeros(npmass,3);
dpm(1,:) = 0;
dpm(2,:) =0;
% ....

% set up linear constraints for clamped wing root
% Number of Degree of freedom
ndof = 3*nnodes;
B = eye(3,ndof);


% retrieve system matrices

[M,K,Z,Qip,f,CRv,CRd,s] = labwing(B, l, b, t, ba, mhinge, rhop, E, G, nelem, dpm);
fprintf("Offset s = %.2f m \n",s);

% compute divergence speed
[udiv,zdiv] = divergence(K, Qip);
fprintf(1,'Divergence speed: %.2f m/s \n', udiv);
return;

