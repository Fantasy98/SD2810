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
mhinge = 3*28e-3;
t = 6E-3;
rhop = 1760;
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
% ....

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

mtot = e1' * M * e1;

# Speed never reach
u_max = 350;
# Stall speed 
u_stall = 70;

S = l * 2*b;
g  = 9.8;
nz = 7;
L = nz * g * mtot;
i = 1;
for u = u_stall:1:u_max
    q = 0.5*Qip.rho * u/3.6* u/3.6;
    CL(i) = L/(q*S);
    i=i+1;
end

loc = CL <= 1.5 ;
loc_u = find(loc==1)(1);

uvec = u_stall:u_max;

umin = uvec(loc_u);



figure(1)
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
x0=15;y0=15;width=600;height=400;


plot([u_stall:u_max] ,CL,"b-","linewidth",3);
hold on 
plot(uvec(loc),CL(loc),"r-","linewidth",3)

plot([u_stall-10  u_max+10],[1.5 1.5],"k-.","linewidth",3)

plot([umin umin],[CL(loc_u) CL(loc_u)],"rp","markersize",20.5,'MarkerFaceColor','red')
leg = legend({
        "Excluded Speed",...
        "Valid Speed",...
        "CL = 1.5",...
        "Minimum Speed"
        });
set(leg,"fontsize",1,"location","northeast");
xlb = xlabel("u (m/s)");
ylb = ylabel("CL");
set([xlb,ylb],"fontsize",20);
% axis([14 24]);
a = get(gca,'XTickLabel');
b = get(gca,'YTickLabel');
set(gca,'XTickLabel',a,'fontsize',20)
set(gca,'YTickLabel',b,'fontsize',20)
 print -djpg umin7.jpg
fprintf("The minumum speed at load factor %.2f is %.2f km/h \n",nz,umin);