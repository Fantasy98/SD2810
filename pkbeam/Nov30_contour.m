
% Nov30 Contour of how mass influence  flutter speed
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

% Here the mhinge is just mass of one rod 
% In the function labwing, it will be tripped
mhinge = 28e-3;

% mhinge = 0;
t = 0.004;%m
rhop = 1950; % Measured Density

% Measured E and G by viberation test
E = 31.5E9;
G = 5.52E9;
ndof = 3*nnodes;
B = eye(3,ndof);

npmass = 8 ;
m1 = (40.33+6.39+2)/1000; %kg
m2 = (20.06+4+2*2)/1000; %kg
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

Nx = 10;
Ny = 10;

x_c = linspace(-b,b-ba,Nx+1);
y_c = linspace(0,l,Ny+1);

x_rev = x_c(end:-1:1);
y_rev = y_c(end:-1:1);
[x_grid,y_grid] = meshgrid(x_c,y_c);
nmode = 3;
neig = 2;
U = zeros(size(y_grid));
     
ex_mass3 = 150/1000;%kg
ex_mass2 = 100/1000;%kg
ex_mass1 = 50/1000;%kg
for i = 1:1:Ny+1
    for j = i:1:Nx+1
    fprintf("At x= %.3f m,y=%.3f m \n",x_c(i),y_c(j))
    dpm(8,:) = [2*ex_mass3 x_c(i) y_c(j)];
    
    [M,K,Z,Qip,f,CRv,CRd,s] = labwing_verbose(B, l, b, t, ba, mhinge, rhop, E, G, nelem, dpm);
    dQip = read_dlm(Z);    
    [Km,Mm,Zm,mQip]=ReduceDim(M,K,dQip,nmode);
    [ucrit,pcrit,zcrit,pconv] = flutter(Mm,Km,mQip,neig);
    U(i,j) = ucrit;
    end
end

% Plot Contour
figure(15)
contourf(x_grid,y_grid,U);
colorbar();
xlab = xlabel("X Coordinate (m)");
ylab = ylabel("Y Coordinate (m)");
set([xlab,ylab],"fontsize",18);
print -djpg DLM_Contour_Flutter_450.jpg
return; 
% retrieve system matrices

% [M,K,Z,Qip,f,CRv,CRd,s] = labwing(B, l, b, t, ba, mhinge, rhop, E, G, nelem, dpm);
% fprintf("Offset s = %.2f m \n",s);

%%% Implement reduced space by introduce only first n modes



fprintf("\nFlutter Speed  is %.2f m/s\n",ucrit);
fprintf("\nFlutter Frequency  is %.2f rad/s\n",pcrit);
% Recover the flutter mode into FEM dimension
z_crit = Zm*zcrit;
% Visualize the mode and rootlocus
vismode(z_crit);
Rootlocus(pconv,neig);