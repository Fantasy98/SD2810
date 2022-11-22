
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

rhop = 1963.7; % Measured Density
fprintf("Reference Mass is 494g \n");

E = 31.5E9;
% Assumed Possion Ratio
G = 5.52E9;


% definition matrix for discrete point masses to attach

% the dpm should give it not only the mass but also the coordinate
npmass = 0 ;
dpm = zeros(npmass,3);
dpm(1,:) = 0;
dpm(2,:) =0;
% ....

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

figure(1)
plot(uvec,deform_tip,"-o","linewidth",1.5);
hold on 
plot([udiv udiv],[min(deform_tip)-10 max(deform_tip)+10])
xlabel("speed (m/s)","fontsize",12);
ylabel("Deflection (m)","fontsize",12);

% for iu = 1:5
%     figure(1+iu)
%     plotmode(V(:,iu));
%     hold on
% end
% legend();
% xlabel("coordinate");
% ylabel("deflection (m)");


return;

