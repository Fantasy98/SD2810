
## Dec14 Final  
% Wing of Full span model under quasi-steady state
% 
% (c) 2004-2016 Dan Borglund <dodde@kth.se> and David Eller <dlr@kth.se>

clear all;

% setup geometry and structural properties
% number of finite elements requested should be a multiple of 3
% since we are simulating the 2 wings, thus the number of element should be 3n
%  The aircraft is like:   ---o---
nelem =24;
nnodes = nelem + 1;

% lab wing dimensions and properties
# Semi-span
l =12.1; % m
% S = 7.41/2; % m^2
c = 0.64;
# According to the technic report, the main wing is assumed as ellipse with chord = 500mm
b = c/2 ; % m
# The alieron span is defined as c*E, E = 0.225 in W12C-JARA-0001.pdf P
ba = c*0.225/2; % m
# mass of hinge need to be figured out 2022/12/12
# Maybe use the same value as labwing is ok
mhinge = 5*28e-3;
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

npmass = 7 ;
m1 = 3*(40.33+6.39+2)/1000; %kg
m2 = 3*(20.06+2*2)/1000; %kg
m3 = 3*(40.33+2*6.39+2*2)/1000; %kg

m_fuselage = 260 - 88;
m_pilot = 53; 
m_parachute = 5;
m_fuse = m_fuselage + m_pilot + m_parachute;
x_coord = (b-ba);
dpm = zeros(npmass,3);
dpm(1,:) = [m1 x_coord 0];
dpm(2,:) =[m2 x_coord 0.5*l/6];
dpm(3,:) =[m3 x_coord  2*0.5*l/6];
dpm(4,:) =[m2 x_coord  3*0.5*l/6];
dpm(5,:) =[m3 x_coord  4*0.5*l/6];
dpm(6,:) =[m2 x_coord 5*0.5*l/6];
dpm(7,:) =[m1 x_coord 0.5*l];

dpm(8,:) = [m_fuse 0 l/2];


dpm(9,:)  =[m2 x_coord  (l/12 + l/2)];
dpm(10,:) =[m3 x_coord  (2*l/12 + l/2)];
dpm(11,:) =[m2 x_coord  (3*l/12 +l/2)];
dpm(12,:) =[m3 x_coord  (4*l/12 + l/2)];
dpm(13,:) =[m2 x_coord  (5*l/12 + l/2)];
dpm(14,:) =[m1 x_coord  l];


ndof = 3*nnodes;
B = [];

% retrieve system matrices
[M,K,Z,Qip,f,CRv,CRd] = nwing(B, l, b, t, ba, mhinge, rhop, E, G, nelem, dpm);




Q0 = Qip.Qtab(:,:,1);


e1 = zeros(ndof,1);
e1(1:3:end) = 1;
e3 = zeros(ndof,1);
e3(3:3:end) = 1;

yn = 0.5*l*linspace(-1,1,nnodes);
B = zeros(3,ndof);
B(1,:) = e1';
B(2,2:3:end) = 1;
B(2,1:3:end) = yn;
% The 
idof = 3*(nelem/2)+3;
B(3,idof) = 1;
Z = null(B);
% Then we derive total mass
mtot = e1' * M * e1;

% Now assume a speed which is lower than the divergence speed
u = 250/3.6;
q = 0.5*Qip.rho*u*u;

g = 9.81;
# Compute the Lift when manuver since L = mg
L = mtot*g;

% Compute Wing Span area
S = l * b * 2;


######### Caution !!! ##############
# For computing flutter, Cla = 2*pi is enough 
# k reduced frequency, for the tail, 
# your k should be adjusted according to chord length of horizontal tail
# assume we are at steady state now, k = 0 
et1 = [1 0 0 1 0 0]';
et3 = [0 0 1 0 0 1]';
k = 0;
# Basic dimension
Stail = 0.86;# Area: fixed
le = 2.15; # tail span fixed, we need full span length as element length in function
bt = 0.5*Stail/le; # semi-chord 
xa = 0.05;
Qetail = beam_amatrix(k, le, bt, xa);
CLtaila = et1' * Qetail * et3/(Stail);

Cmtaila = et3' * Qetail * et3/(Stail*2*bt);

xt = 2.67;
T2 = [  1 -le/2 -xt
        0    1   0  
        0    0   1
        1  le/2  -xt 
        0    1   0 
        0    0   1];

# Make matrix has same size as Q0
Qtail = zeros(size(Q0));

#Find location of fuselarge nodes
# 
ivec = [3*fix(nnodes/2)+1 : 3*fix(nnodes/2)+3];

# Give the contribution of tail 
# Where T1 could be T2'
Qtail(ivec,ivec) = T2' * Qetail * T2;
# consider the stablizer now  


Q = Q0 + Qtail;

% Qa = zeros(size(Q0));
% for i = 1:ndof
%         Qa(:,i) = f';
        
% end


% Q1 = Q0 + Qtail + Qa;




% Qaz = zeros(size(Q0));
% for i = 1:ndof
%         Qaz(:,i) = CRv;
% end
% Q00 = Q0 + Qtail + Qaz';

% Qa1 = zeros(size(Q0));
% for i = 1:ndof
%         Qa1(:,i) = f'.*CRd;
% end
% Q11 = Q0 + Qtail + Qa1';


CLa = e1' * Q0 * e3/(S);
Cma = e3' * Q0 * e3/(S*2*b);

CLde = e1' * Qtail * e3/(S);
Cmde = e3' * Qtail * e3/(S*2*b);

############# Ineria Relief
% Make elastic deformation is orthorgonal to the rigid body
% Linear function distribution of right slope

% Check the rank 
ndof == rank([B' Z]);


u = 250/3.6;
q = 0.5*Qip.rho*u*u;

####Solve for initial angle at inerial relief condition,
# consider pitching moment equilibrium
 
        
% LHS =[  Z'*(K-q*Q1)*Z   -q*Z'*Q*e3    -q*Z'*Q1*e3

%         q*e1'*Q1*Z      q*e1'*Q*e3    q*e1'*Q1*e3
%         % q*e1'*Q1*Z      q*e1'*Q*e3    q*Z'*f

% %        0.07*q*e3'*Q*Z     0.07*q*e3'*Q*e3    (ba+b+0.07)*q*e3'*f ];
%         % (0.5*b-0.07)*q*e3'*Q*Z     (0.5*b-0.07)*q*e3'*Q*e3    (0.5*b-(ba+b+0.07))*q*e3'*f ];
%         q*e3'*Q1*Z/(S*c)         q*e3'*Q*e3/(S*c)         q*e3'*Q1*e3/(S*c) ];
%         % q*e3'*Q1*Z         q*e3'*Q1*e3         CRd ];
%         % q*e3'*Q*Z         q*e3'*Q*e3         -CRd ];

LHS = [
        Z'*(K-q*Q)*Z    -q*Z'*Q*e3      -q*Z'*f

        q*e1'*Q*Z       q*e1'*Q*e3      q*e1'*f

        q*e3'*Q0*Z       q*e3'*Q0*e3      q*e3'*f
        ];


nz = 7;

RHS = [-nz*g*Z'*M*e1

        nz*g*e1'*M*e1

        0
        ];

x = LHS\RHS;



alfa0  = x(end-1);
delta0 = x(end);
vhat   = x(1:end-2);
vtot =  Z*vhat + alfa0*e3 + delta0*e3;
va = Z*vhat + alfa0*e3;
res = K*vtot + M*nz*g*e1 - q*Q0*vtot;

fprintf("Deformation computed, the residual =  %E \n",min(abs(res)))
rank = rank([B' Z]);

alfa = vtot(end)-delta0;
alfae = vtot(end)-delta0-alfa0;

fprintf("At speed = %.2f m/s \n",u);
fprintf("The angle of attack now is = %.2f deg \n",alfa*180/pi)
fprintf("The rigid angle of attack  = %.2f deg \n",alfa0*180/pi)
fprintf("The elastic angle of attack  = %.2f deg \n",alfae*180/pi)
fprintf("The aileron  is = %.2f deg \n",delta0*180/pi)
% fprintf("CL is = %.2f\n",CL)
# Plot the stress distribution 
plot_stress(Z*vhat,l,b,t,E,G,c);
