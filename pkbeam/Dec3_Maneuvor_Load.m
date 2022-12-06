
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
% Then we derive total mass
mtot = e1' * M * e1;
% Now assume a speed which is lower than the divergence speed
u = 15;
q = 0.5*Qip.rho*u*u;

g = 9.81;
# Compute the Lift when manuver since L = mg
L = mtot*g;


% For squared wing, Lift coefficient could be assumed as 2pi
CLalfa  = 2*pi;
% Compute Wing Span area
S = l * b * 2;
% To compute angle of attack
% Assume the rigid body case, we do not take elastic defromation of wing into consideration
% Lift = q*S*Cla*alfa and Lift = Mass

alfa = L/(q*S*CLalfa);
alfa = alfa*180/pi;
fprintf("Lift rigid body only = %.2f N \n",L);
fprintf("AOA rigid body only = %.2f deg \n",alfa);

# Try to match the lift by adjusting AOA
% v0 represents reference elastic deformation alpha0, which is at 3rd dof
e3 = zeros(ndof,1);
e3(3:3:end) = 1;
alfa1 = 0.096;
v0 = alfa1 * e3;
% By solving the governing equation to get deformation in null space 
vehat = (Z'*(K-q*Q0)*Z)\(q*Z'*Q0*v0);
% Then the total deformation is the combine of v0 Z*vehat (recover from null space)
vtot =v0 + Z*vehat ;
Lflex = q*e1'* Q0 * vtot;
fprintf("Now if we have AOA = %.2f deg \n",alfa1*180/pi);
fprintf("Considering elastic deformation, Lift is %.2f N \n",Lflex);



#############################################
% To derive the elastic deformation and Intial angle of attack 
% Now assume we are at steady state, where No accleration for M 
% 1. From Governing equation in null space we get 
%       Z'* (K - q * Q0) * Z * vehat - q* Z'*Q0 *e3.*alpha0 = 0 
% 2. From quasi-steady flight condition, we have Lift = Weight
%       q*e1'*Q0 * Z * vehat + q*e1'*Q0*e3.*alpha0 = mg 
% Combing 1&2 we close the equation and written into matrix form:

% Left hand side
LHS = [ Z' *(K-q*Q0)*Z  -q*Z'*Q0*e3
        q*e1'*Q0*Z       q*e1'*Q0*e3];

% Right hand side
RHS = [zeros(ndof-3,1) 
        mtot*g];

% Solving [vehat
%          alpha0]
x = LHS\RHS;
% Then we now the last element in the vector is initial angle of attack

alfa0 = x(end)*180/pi;
fprintf("If Elastic deformation, the initial angle of attack is %.2f deg \n",alfa0);
# Compare with the alfa0 from rigid body problem
% IF we increase stiffness , the angle of attack will be equal to the rigid body case
K_enhance = 10000*K;
% Left hand side
LHS = [ Z' *(K_enhance-q*Q0)*Z  -q*Z'*Q0*e3
        q*e1'*Q0*Z       q*e1'*Q0*e3];

% Right hand side
RHS = [zeros(ndof-3,1) 
        mtot*g];
% Retrieve [vehat
%          alpha0]
x = LHS\RHS;
% Then we now the last element in the vector is initial angle of attack
alfa0 = x(end)*180/pi;
fprintf("If elastic deformation with enhance stiffness ,alfa0 =  %.2f deg \n",alfa0);


# Visualize the mode by plotmode
# Now you can see the initial twist angle 
alfa0 = x(end);
vehat = x(1:end-1);
vtot = alfa0* e3 + Z*vehat;
plotmode(vtot)

###########################################
# Inertia relief 
# New B matrix
B = eye(3,ndof);
B(1,1:3:end) = 1;
Z = null(B);
vacc = g*e1;
%K*v+M*vacc = q*Q0*v
norm(K*B(1,:)')



return;
