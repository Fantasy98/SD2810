
% Dec 1 
% Improve Flutter based on DLM method 
% Propose at least a way to imporve the flutter speed 
% Not exceed 600g !!!
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

nmode = 6;
neig = 3;

npmass = 12 ;
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
########### Extra Mass ###########
ex_mass1 = 50/1000;
ex_mass2 = 100/1000;
ex_mass3 = 150/1000;
% Method 1 Flutter speed = 24.4/s for DLM & 20.50 for Strip

dpm(8,:) = [2*ex_mass3,-b l];
dpm(9,:) = [2*ex_mass2,-b (5.5/6) *l]; % 1.4667m

% % Method 2 Flutter speed = 24.9m/s for DLM & 20.90m/s for Strip
dpm(8,:) = [ex_mass3+ex_mass3,-b l];
dpm(9,:) = [ex_mass2+ex_mass2,-b (4.8/5) *l]; %1.5360m

% % Method 3 Flutter speed = 25.50 m/s for DLM & 21.4m/s for Strip
dpm(8,:) = [ex_mass3+ex_mass3,-b l];
dpm(9,:) = [ex_mass2+ex_mass2,-b (4.8/5) *l];
dpm(10,:) = [ex_mass1*2,-b (3.5/5) *l]; % 1.12m

% % Method 4 Flutter speed = 25.700 m/s for DLM & 21.50m/s for Strip
dpm(8,:) = [ex_mass3+ex_mass3,-b l];
dpm(9,:) = [ex_mass2+ex_mass2,-b (4.8/5) *l];
dpm(10,:) = [ex_mass1*2,-b (3.7/5) *l]; % 1.184m


##################################

[M,K,Z,Qip,f,CRv,CRd,s] = labwing_verbose(B, l, b, t, ba, mhinge, rhop, E, G, nelem, dpm);
dQip = read_dlm(Z);
[Km,Mm,Zm,mQip]=ReduceDim(M,K,dQip,nmode);
[ucrit,pcrit,zcrit,pconv,uvec] = flutter(Mm,Km,mQip,neig,150);
fprintf("\nFlutter Speed  is %.3f m/s\n",ucrit);
% fprintf("Flutter Frequency  is %.3f rad/s\n",pcrit);
[udiv,zdiv] = divergence(K, dQip);
fprintf(1,'Divergence speed: %.3f m/s \n', udiv);
[urev,zrev] = reversal(K, dQip, f, CRv, CRd);
% fprintf(1,'Reversal speed: %.3f m/s \n', urev);
% Rootlocus(pconv,neig);
% figure(14)
plot(uvec,real(pconv),"o");



#######################################################################################################################
fprintf("\n Strip Theory Results \n")
% Strip Theory
[M,K,Z,Qip,f,CRv,CRd,s] = labwing_verbose(B, l, b, t, ba, mhinge, rhop, E, G, nelem, dpm);

[Km,Mm,Zm,mQip]=ReduceDim(M,K,Qip,nmode);

[ucrit,pcrit,zcrit,pconv,uvec] = flutter(Mm,Km,mQip,neig);
fprintf("\nFlutter Speed  is %.2f m/s\n",ucrit);
fprintf("Flutter Frequency  is %.2f rad/s\n",pcrit);
[udiv,zdiv] = divergence(K, Qip);
fprintf(1,'Divergence speed: %.2f m/s \n', udiv);
%compute reversal speed
[urev,zrev] = reversal(K, Qip, f, CRv, CRd);
fprintf(1,'Reversal speed: %.2f m/s \n', urev);
return



figure(5)
Rootlocus(pconv,neig);
% print -djpg Strip_rootlocus.jpg




figure(12)
for imode = 1:3
    plot(uvec, real(pconv(imode,:)),"o-","linewidth",0.8,"markersize",3.5);
    hold on 
    end
plot([10 25],[0,0],"k-.","linewidth",1.5);
leg = legend({
        "Real part of Mode 1",...
        "Real part of Mode 2",...
        "Real part of Mode 3",...
        "Re(p) = 0"
        });
set(leg,"fontsize",8,"location","southeast");
xlabel("u (m/s)");
ylabel("Real(p)");
axis([14 25]);
% print -djpg Strip_Rep_u.jpg

% % Recover the flutter mode into FEM dimension
% z_crit = Zm*zcrit;
% % Visualize the mode and rootlocus
% vismode(z_crit);
% Rootlocus(pconv,neig);


##############################################################################################################
