
% For test the frequency with analytical eigenfrequncy
% Free End Beam
clear all;

% setup geometry and structural properties
% number of finite elements requested should be a multiple of 3
nelem = 24;

nnodes = nelem + 1;

% lab wing dimensions and properties
l = 1.63; % m
b = 0.28/2;
% Here we do not consider the aileron 
%and mass of hinge since the test will be apply on the flat plate
ba = 0;
mhinge = 0;
t = 0.004;
% The estimate data from David
% which should be adjust after experiment !
%################################
mass = 3.585; % kg
volume =l * 2*b * t;
rhop = mass/volume;
fprintf("Test density is %.2f \n",rhop);


% A guess of E and G 
E = 25E9;
E = 31.5E9;
% Assumed Possion Ratio
% E = E * 2.5;
possion = 0.21;

G = E/2*(1+possion);
G = 5.52E9;

%################################
% From Dave
% E = (25.944+5.4)*10^9;
% G = 5.53E9;
% definition matrix for discrete point masses to attach
npmass = 0 ;
dpm = zeros(npmass,3);
dpm(1,:) = 0;
dpm(2,:) =0;
% ....

% set up linear constraints for clamped wing root
% Number of Degree of freedom
ndof = 3*nnodes;
% Clamped
% B = eye(3,ndof);
% No constrain == Free
B = [];
 
% Here we are testing condition without clamping

% retrieve system matrices
[M,K,Z,Qip,f,CRv,CRd,s] = labwing(B, l, b, t, ba, mhinge, rhop, E, G, nelem, dpm);

fprintf("Offset s = %.2f m \n",s);

% Compute the Inneria
I= (2*b*t^3)/12;

% print free vibration frequencies
[V,LAMBDA] = eig(K,M);
Vhat = Z * V;
omega = diag(sqrt(LAMBDA)./(2*pi));
freq = real(omega(1:20));
for i = 1:4
    fprintf(" \n Frequency-%.2d:  \n %.2frad/s \n %.2frad/s \n %.2frad/s \n",i,freq(1+3*i),...
                                                                                    freq(2+3*i),...
                                                                                    freq(3+3*i))

    % fprintf(" 2 + 3N Bending frequencies is %.2f rad/s \n",freq(2+3*i))

    % fprintf(" 3 + 3N Torsional frequencies is %.2f rad/s \n",freq(3+3*i))
end

fprintf("\n E=%.3E \n",E);
fprintf("\n G=%.3E \n",G);
fprintf("\n Material Density rho=%.4d kg/m^3 \n",rhop);

% Theoritical Formula for frequncies
% fprintf("Analytical Solution of Frequencies: \n")
% 1. Bending frequency: 
% Compute beam mass per length
my = rhop * 2*b*t;
% Constant Bn values for  3 lowest modes for free
% uL = [0 4.730 7.853];
% r = 3:10;
% url = (2.*r+1/2).*pi;
% UL = [uL url];
% f_bending = ((UL).^2)*sqrt(E*I/(my*l^4));
% %fprintf("Frequency of bending is %.2f rad/s \n", f_bending );

return;

