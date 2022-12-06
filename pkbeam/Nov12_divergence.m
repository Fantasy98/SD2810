
% Programme for computing divergence speed
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
% ba = 0;
% measured from lab
% mhinge = ( 2*(40.33+6.39+2)+...
%             2*(20.06+2*2)+...
%             2*(40.33+2*6.39+2*2)+...
%             3*28.9)  % g
% mhinge = mhinge/1000; % kg 
%% mhinge is only one rode
mhinge = 28e-3;
t = 0.004;%m

rhop = 1950; % Measured Density


E = 31.5E9;
% Assumed Possion Ratio
G = 5.52E9;


% definition matrix for discrete point masses to attach

% the dpm should give it not only the mass but also the coordinate
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

% set up linear constraints for clamped wing root
% Number of Degree of freedom
ndof = 3*nnodes;
B = eye(3,ndof);


% retrieve system matrices

[M,K,Z,Qip,f,CRv,CRd,s] = labwing(B, l, b, t, ba, mhinge, rhop, E, G, nelem, dpm);
fprintf("Offset s = %.2f m \n",s);

% Derive the deflection VS speed
% set angle of attack a0 = 1 
% Steady State aerodynamic force
%[K-qQ]v = qQv0 where v0 has the initial angle of attack
% Steady aerodynamic force
A = Qip.Qtab(:,:,1);
% angle of attack
vshape = size(K,1);
% Initial angle of attack alpha0 = 1 deg
v0 = zeros(vshape,1);
v0(3:3:end) = pi/180;
iu = 0;
for u = 19:0.1:21
    iu = iu+1;
    % dynamic pressure
    q = 0.5*Qip.rho * u *u
    % solve the deflection 
    vt = (K - q*A)\(q*A*v0);
    v = vt + v0;
    % store all deflection
    V(:,iu) = v;
    deform_tip(iu) = v(end-2);
    uvec(iu) = u;
end

% compute divergence speed
[udiv,zdiv] = divergence(K, Qip);
fprintf(1,'Divergence speed: %.2f m/s \n', udiv);

[maxp,loc] = max(deform_tip);
figure(1)
plot(uvec,deform_tip,"-o","linewidth",3.5,"markersize",12.5);
hold on 
plot([uvec(loc) uvec(loc)],[min(deform_tip)-40 max(deform_tip)+40],"-.","linewidth",4.5)
axis([min(uvec) max(uvec) min(deform_tip)-15 max(deform_tip)+15  ]);
leg = legend({"deformation at wing tip","Divergence Speed"})
set(leg,"location","northwest","fontsize",20);
xlab=xlabel("speed (m/s)","fontsize",13);
ylab=ylabel("Deformation (m)","fontsize",13);
set([xlab ylab],"fontsize",20)
a = get(gca,'XTickLabel');
b = get(gca,'YTickLabel');
set(gca,'XTickLabel',a,'fontsize',20)
set(gca,'YTickLabel',b,'fontsize',20)
print -djpg -r0 Divergence_vs_speed.jpg
% for iu = 1:5
%     figure(1+iu)
%     plotmode(V(:,iu));
%     hold on
% end
% legend();
% xlabel("coordinate");
% ylabel("deflection (m)");


return;

