
% Nov 14 : Aileron Reversal
clear all;

% setup geometry and structural properties
% number of finite elements requested should be a multiple of 3
nelem = 24;
nnodes = nelem + 1;

% lab wing dimensions and properties
l =1.6; % m
b = 0.175; % m
ba = 0.03 ;
% measured from lab % Not correct need to be calculate again
mhinge = 0;
t = 0.004;%m

rhop = 1950; % Measured Density

% Measured E and G by viberation test
E = 31.5E9;
% Assumed Possion Ratio
G = 5.52E9;


% definition matrix for discrete point masses to attach
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

% set deflection of alieron 
% The certia, set 
delta = 5*pi/180; % rad/s
% Set speed
iu = 0; 
for u = 15:.1:20
iu = iu+1;
q = 0.5* Qip.rho * u * u;
A = Qip.Qtab(:,:,1);

% Solve for the deformation
v = (K-q*A)\(q*f*delta);
    
    % plotmode(v);

    % Plot deformation 
v_all(iu,:) = v;
deform_tip(iu) = v(end-2);
uv(iu) = u;
end
% Plot the deformation VS velocity 
figure(1)
plot(uv,deform_tip,"linewidth",1.5);
xlabel("Speed (m/s)","fontsize",12);
ylabel("Deformation (m)","fontsize",12);


% Compute the aileron reverse
A_eig = A - f*CRv./CRd;
[V,D] = eig(A_eig , K);
qrev = 1/max(diag(D));
urev = sqrt(2*qrev/Qip.rho);

figure(2)
set(get(gca, 'Title'), 'String', 'The deformation before reversal speed');
plotmode(v_all(30,:));

figure(3)
set(get(gca, 'Title'), 'String', 'The deformation after reversal speed');
plotmode(v_all(42,:));

% The reversal speed should be less than divergence
% The flutter velocity should be lower than the divergence
% u_flutter = 16 
% Flutter velocity is expected to be lower than reversal in the begining
% Then we can improve it by adding concentrated mass on the wing
% To put mass and imporve the flutter speed.
[urev,zrev] = reversal(K,Qip,f,CRv,CRd);
fprintf("The reverse velocity is %.2f\n",urev);
figure(4)
set(get(gca, 'Title'), 'String', 'The deformation at reversal speed');

plotmode(zrev);
