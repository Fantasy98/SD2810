clear all
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

% definition matrix for discrete point masses to attach
% All attachment of contorl surface should be assigned as dpm
% The format is dpm(i,:) = [mass x_coord y_coord]
% 7 connection in total
% From  root to tip = 1:7
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

dpm(8,:) =[m_fuse 0 l/2];


dpm(9,:)  =[m2 x_coord  (l/12 + l/2)];
dpm(10,:) =[m3 x_coord  (2*l/12 + l/2)];
dpm(11,:) =[m2 x_coord  (3*l/12 +l/2)];
dpm(12,:) =[m3 x_coord  (4*l/12 + l/2)];
dpm(13,:) =[m2 x_coord  (5*l/12 + l/2)];
dpm(14,:) =[m1 x_coord  l];

ndof = 3*nnodes;
% Since aircraft is maneuvoring, no rigid body conrtain at all
% No constrain before we get M, K ,Qip 
B = []

% retrieve system matrices
[M,K,Z,Qip,f,CRv,CRd] = nwing(B, l, b, t, ba, mhinge, rhop, E, G, nelem, dpm);


e1 = zeros(ndof,1);
e1(1:3:end) = 1;
e3 = zeros(ndof,1);
e3(3:3:end) = 1;

## Assemble constrain matrix
yn = 0.5*l*linspace(-1,1,nnodes);
B = zeros(3,ndof);
B(1,:) = e1';
B(2,2:3:end) = 1;
B(2,1:3:end) = yn;
idof = 3*(nelem/2)+3;
B(3,idof) = 1;

nz=5.3;
k = 0;
i = 1;
for ui = 172:5:350
    u = ui/3.6;
    [Q,vtot,ve,alfa0,delta0] = flight_load(K,M,B,Qip,f,nelem,u,nz,k);
    [SF_n,SF_s] = safety_factor(vtot);
    SFN5(i) = SF_n;
    SFS5(i) = SF_s;
    i = i+1;
end


nz=7;
k = 0;
i = 1;
for ui = 172:5:350
    u = ui/3.6;
    [Q,vtot,ve,alfa0,delta0] = flight_load(K,M,B,Qip,f,nelem,u,nz,k);
    [SF_n,SF_s] = safety_factor(vtot);
    SFN7(i) = SF_n;
    SFS7(i) = SF_s;
    i = i+1;
end

uall = 172:5:350;
figure(1)
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
x0=15;y0=75;width=600;height=400;
set(gcf,'units','points','position',[x0,y0,width,height])
plot(uall,SFN5,"-ob","linewidth",1.8);
hold on 
plot(uall,SFN7,"-og","linewidth",1.8);
plot(uall,SFS5,"-or","linewidth",1.8);
plot(uall,SFS7,"-oy","linewidth",1.8);

leg = legend({
                "Normal Stress at +5.3"
                "Normal Stress at +7"
                "Shear Stress at +5.3"
                "Shear Stress at +7"

                    });
set(leg,"fontsize",15,"location","northwest")

ylab =ylabel("Safety factor","fontsize",8)
xlab =xlabel("Speed (km/h)","fontsize",8)
set([xlab,ylab],"fontsize",15)
a = get(gca,'XTickLabel');
b = get(gca,'YTickLabel');
set(gca,'XTickLabel',a,'fontsize',15)
set(gca,'YTickLabel',b,'fontsize',15)

print -djpg SFfactor.jpg