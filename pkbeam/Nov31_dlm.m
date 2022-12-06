
% Nov31 
% load dlm data for 24 element 
% Improve 
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

x_coord_extra = -(x_coord+0.05);
m_extra = [50 100 150]/1000;

%% A extreme case which will trigger the mode3 be hump as well 
%% In this way, we are decreasing flutter speed instead of improve it !
% dpm(9,:) = [3*m_extra(3),-b,l];



[M,K,Z,Qip,f,CRv,CRd,s] = labwing(B, l, b, t, ba, mhinge, rhop, E, G, nelem, dpm);

%Dlm 
dQip = read_dlm(Z);
nmode = 6;
neig = 3;
[Km,Mm,Zm,mQip]=ReduceDim(M,K,dQip,nmode);

[ucrit,pcrit,zcrit,pconv,uvec] = flutter(Mm,Km,mQip,neig);

fprintf("\nFlutter Speed  is %.2f m/s\n",ucrit);
fprintf("\nFlutter Frequency  is %.2f rad/s\n",pcrit);

% % compute divergence speed
[udiv,zdiv] = divergence(K, dQip);
fprintf(1,'Divergence speed: %.2f m/s \n', udiv);
%compute reversal speed
[urev,zrev] = reversal(K, dQip, f, CRv, CRd);
fprintf(1,'Reversal speed: %.2f m/s \n', urev);
return
figure(4)
Rootlocus(pconv,neig);
print -djpg DLM_rootlocus.jpg

figure(11);
for imode = 1:3
    plot(uvec, real(pconv(imode,:)),"o-","linewidth",0.8,"markersize",4.5);
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
print -djpg DLM_Rep_u.jpg




% Strip Theory
[Km,Mm,Zm,mQip]=ReduceDim(M,K,Qip,nmode);

[ucrit,pcrit,zcrit,pconv,uvec] = flutter(Mm,Km,mQip,neig);
fprintf("\nFlutter Speed  is %.2f m/s\n",ucrit);
fprintf("\nFlutter Frequency  is %.2f rad/s\n",pcrit);
[udiv,zdiv] = divergence(K, Qip);
fprintf(1,'Divergence speed: %.2f m/s \n', udiv);
%compute reversal speed
[urev,zrev] = reversal(K, Qip, f, CRv, CRd);
fprintf(1,'Reversal speed: %.2f m/s \n', urev);
figure(5)
Rootlocus(pconv,neig);
print -djpg Strip_rootlocus.jpg

figure(12)
for imode = 1:3
    plot(uvec, real(pconv(imode,:)),"o-","linewidth",0.8,"markersize",4.5);
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
print -djpg Strip_Rep_u.jpg

return;
% % Recover the flutter mode into FEM dimension
% z_crit = Zm*zcrit;
% % Visualize the mode and rootlocus
% vismode(z_crit);
% Rootlocus(pconv,neig);
