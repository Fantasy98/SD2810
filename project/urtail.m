
## Dec14 Final  
% Wing of Full span model under quasi-steady state
% 
% (c) 2004-2016 Dan Borglund <dodde@kth.se> and David Eller <dlr@kth.se>

clear all;

% setup geometry and structural properties
% number of finite elements requested should be a multiple of 3
% since we are simulating the 2 wings, thus the number of element should be 3n
%  The aircraft is like:   ---o---
nelem =18;
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
x_coord = (280-175)/1000;
dpm = zeros(npmass,3);
dpm(1,:) = [m1 x_coord 0];
dpm(2,:) =[m2 x_coord 27/100];
dpm(3,:) =[m3 x_coord  53/100];
dpm(4,:) =[m2 x_coord  80/100];
dpm(5,:) =[m3 x_coord  106/100];
dpm(6,:) =[m2 x_coord 133/100];
dpm(7,:) =[m1 x_coord 160/100];
dpm(8,:) = [280 0 l/2]
% ....

ndof = 3*nnodes;
% Since aircraft is maneuvoring, no rigid body conrtain at all
% No constrain before we get M, K ,Qip 
B = []

% retrieve system matrices
[M,K,Z,Qip,f,CRv,CRd] = nwing(B, l, b, t, ba, mhinge, rhop, E, G, nelem, dpm);


% First let's compute divergence speed as reference.
% Divergence speed should be identicall to any case since it's not related to the accelration 
####################################
% Obtain static aerodynamic forces
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
u = 350/3.6;
q = 0.5*Qip.rho*u*u;

g = 9.81;
# Compute the Lift when manuver since L = mg
L = mtot*g;

% Compute Wing Span area
S = l * b * 2;
S = 7.41;
% should be 2pi
CLalfa  = e1' * Q0 * e3/(S);

# Shift the ba proportion of b it will change the location of elastic axis
# for example ba = 0.5*b,  Cma = 0!
Cma = e3'*Q0*e3/(S*2*b);


alfa = L/(q*S*CLalfa);

######### Caution !!! ##############
# For computing flutter, Cla = 2*pi is enough 
# k reduced frequency, for the tail, 
# your k should be adjusted according to chord length of horizontal tail
# assume we are at steady state now, k = 0 
et1 = [1 0 0 1 0 0]';
et3 = [0 0 1 0 0 1]';
k = 0;
# Arbitary decide le and b at this stage
le = 2.644;
bt = 0.17 
# xa is the
    % local position of the elastic axis relative the the
    % midchord (xea - xmid)
xa = 0;
# if bt = 0.2 , xa= -0.1 is at 1/4 c which turns Cmtaila = 0 !!
xa = -0.1;
Qetail = beam_amatrix(k, le, bt, xa);

# Surface
Stail=le*2*bt

#Expected to get 2pi as well
CLtaila = et1' * Qetail * et3/(Stail)

Cmtaila = et3' * Qetail * et3/(Stail*2*bt)

xt = 2;
T2 = [  1 -le/2 -xt
        0    1   0  
        0    0   1
        1  le/2  -xt 
        0    1   0 
        0    0   1]

# Make matrix has same size as Q0
Qtail = zeros(size(Q0));

#Find location of fuselarge nodes
# 
ivec = [3*fix(nnodes/2)+1 : 3*fix(nnodes/2)+3]

# Give the contribution of tail 
Qtail(ivec,ivec) = T2' * Qetail * T2;

CLa = e1' * Qtail * e3/(Stail)

Cma = e3' * Qtail * e3/(Stail*2*b)

# consider the stablizer now  
Q0 = Q0 + Qtail;
CLa = e1' * Q0 * e3/(S)
Cma = e3' * Q0 * e3/(S*2*b)

CLde = e1' * Qtail * e3/(S)
Cmde = e3' * Qtail * e3/(S*2*b)

# A way to compute Roll moment 
R =B(2,:)*Q0 * e3

A = [CLa CLde 
        Cma Cmde];


rhs = [mtot*g/(q*S)
        0];

x = A\rhs


############# Ineria Relief
% Make elastic deformation is orthorgonal to the rigid body
% Linear function distribution of right slope

% Check 
ndof == rank([B' Z])

####Solve for initial angle 
LHS = [ Z' *(K-q*Q0)*Z  -q*Z'*Q0*e3
        q*e1'*Q0*Z       q*e1'*Q0*e3];

% Right hand side
nz = 9;
RHS = [-nz*g*Z'*M*e1 
        g*e1'*M*e1];

% Retrieve [vehat
%          alpha0]
x = LHS\RHS;
alfa0 = x(end)*180/pi;
alfa0 =x(end);
vhat = x(1:end-1);
vtot =  Z*vhat + alfa0*e3;
res = K*vtot + M*nz*g*e1 - q*Q0*vtot
fprintf("Deformation computed, bthe residual of results %E \n",min(abs(res)))
rank = rank([B' Z]);

ev = eig(Z'*Q0*Z,Z'*K*Z);
qdiv = 1/max(ev);
udiv = sqrt(2*qdiv/Qip.rho)
plot_stress(vtot,l,b,t,E,G,c)

% K * vtot  + g*M*e1