
% Dec3 Flight Loads 
% Quasi-Steady
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

% Here the mhinge is just mass of one rod 
% In the function labwing, it will be tripped
mhinge = 28e-3;

% mhinge = 0;
t = 0.004;%m
rhop = 1950; % Measured Density

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

ndof = 3*nnodes;
B = eye(3,ndof);
% Since aircraft is maneuvoring, no rigid body conrtain at all
% No constrain before we get M, K ,Qip 
B = []

% retrieve system matrices
[M,K,Z,Qip,f,CRv,CRd,s] = labwing_verbose(B, l, b, t, ba, mhinge, rhop, E, G, nelem, dpm);
fprintf("Offset s = %.2f m \n",s);

%
% Obtain static aerodynamic forces
Q0 = Qip.Qtab(:,:,1);
% Now we compute constrain maunally and convert into null space
B = eye(3,ndof);
Z = null(B);
% Solving divergence speed in null space where K is equal to K we used to use since labwing gives out the K already in null space
% Same trick: Put non-singular matrix at the second place A-1/lambda K = 0 rather than K - lambda A = 0
ev = eig(Z'*Q0*Z,Z'*K*Z);
% Dynamic pressure
qdiv = 1/max(ev);
% Divergence speed
udiv = sqrt(2*qdiv/Qip.rho);

#########################################
% Angle of attack 
% To derive the total mass 
% e1 is  where the gravity accerlation 
% At quasi-steady case, the rotation acceleration is constant 
% When we assume the load factor =1 , it will be equal to g 
% since the gravity will only affect twist dof which is at 3rd dof
%  ve = g.*[ 0 0 1 0 0 1 ...]

e1 = zeros(ndof,1);
e1(1:3:end) = 1;
e3 = zeros(ndof,1);
e3(3:3:end) = 1;
% Then we derive total mass
mtot = e1' * M * e1;
% Now assume a speed which is lower than the divergence speed
u = 15;
q = 0.5*Qip.rho*u*u;

g = 9.81;
# Compute the Lift when manuver since L = mg
L = mtot*g;
CLalfa  = 2*pi;
% Compute Wing Span area
S = l * b * 2;
% To compute angle of attack
% Assume the rigid body case, we do not take elastic defromation of wing into consideration
% Lift = q*S*Cla*alfa and Lift = Mass
alfa = L/(q*S*CLalfa);

############# Ineria Relief
% Make elastic deformation is orthorgonal to the rigid body
B = zeros(3,ndof);
B(1,:) = e1';
B(2,2) = 1;
B(3,3) = 1;
Z = null(B);
% Check 
ndof == rank([B' Z])

####Solve for initial angle 
LHS = [ Z' *(K-q*Q0)*Z  -q*Z'*Q0*e3
        q*e1'*Q0*Z       q*e1'*Q0*e3];

% Right hand side
nz = 1;
RHS = [-nz*g*Z'*M*e1 
        g*e1'*M*e1];
% Retrieve [vehat
%          alpha0]
x = LHS\RHS;
x(end)*180/pi;

% Divergence
ev = eig(Z'*Q0*Z,Z'*K*Z);
qdiv = max(1/ev);
udiv = sqrt(2*qdiv/Qip.rho);