% Nov 23 Improve the flutter speed by adding discrete mass on the wing
% mhinge has been updated
% pk_bounds and pk_bisect

% Updated derivation of flutter speed
clear all;

% setup geometry and structural properties
% number of finite elements requested should be a multiple of 3
nelem = 15;
nnodes = nelem + 1;
l =1.6; % m
b = 0.175; % m
ba = 0.03 ;

% Here the mhinge is just mass of one rod 
% In the function labwing, it will be tripped
mhinge = 28e-3;
t = 0.004;%m

rhop = 1950; % Measured Density

% Measured E and G by viberation test
E = 31.5E9;
% Assumed Possion Ratio
G = 5.52E9;


% definition matrix for discrete point masses to attach
npmass = 12 ;

m1 = (40.33+6.39+2)/1000; %kg
m2 = (20.06+2*2)/1000; %kg
m3 = (40.33+2*6.39+2*2)/1000; %kg
% According to the photo
x_coord = (280-175)/1000;

% extra mass we can use to improve flutter speed
m_extra = [50 100 150]/1000;

% Fixed mass
dpm = zeros(npmass,3);
dpm(1,:) = [m1 x_coord 0];
dpm(2,:) =[m2 x_coord 27/100];
dpm(3,:) =[m3 x_coord  53/100];
dpm(4,:) =[m2 x_coord  80/100];
dpm(5,:) =[m2 x_coord  106/100];
dpm(6,:) =[m2 x_coord 133/100];
dpm(7,:) =[m2 x_coord 160/100];

% To improve flutter speed
% Adding extra discrete mass on the wing 
x_coord_extra = x_coord + 0.06;
% x_coord_extra = b/2;
dpm(8,:) = [m_extra(2) -x_coord_extra 0];
dpm(9,:) = [m_extra(2) -x_coord_extra l/4];
dpm(10,:) = [m_extra(2) -x_coord_extra l/2];
dpm(11,:) = [m_extra(2) -x_coord_extra 3*l/4];
dpm(12,:) = [m_extra(2) -x_coord_extra l];
% Test before/after adding hingmass 
% In theory, since hinges are after the elastic axis, the flutter speed should be decreased
% dpm = zeros(npmass,3);

% Consider discrete mass of hinges, the flutter speed is 15.7
% Without mass of hinges, the flutter speed is 16.0


% set up linear constraints for clamped wing root
% Number of Degree of freedom
ndof = 3*nnodes;
B = eye(3,ndof);

% retrieve system matrices

[M,K,Z,Qip,f,CRv,CRd,s] = labwing(B, l, b, t, ba, mhinge, rhop, E, G, nelem, dpm);
fprintf("Offset s = %.2f m \n",s);



%#######################################
% 1 Check flutter speed and compare with the divergence and reversal 
% Note that ieig is the index of eigvalue, 
% which means it can only solve for single mode 
u = 15.5; 
i = 0.05;
neig = 3;
for iu = 1:100
    
    [kbounds] = pk_bounds(u,M,K,Qip,neig);    
    for ieig = 1:neig
        [pu,vu] = pk_bisect(ieig,u,M,K,Qip,kbounds);
        pconv(ieig,iu) = pu;
        Vu(ieig,iu,:) = vu;
    end
    
    u = u + i;
    uvec(iu) = u;
end



% Plot Root Locus
figure(1);
for i = 1:neig
plot(real(pconv(i,:)),imag(pconv(i,:)),"o-","linewidth",0.8);
hold on 
end
plot([0 0],[0 1],"k-","linewidth",1.5);
axis([-0.2 0.2 0 0.5]);
leg = legend({
        "Root Locus of Mode 1",...
        "Root Locus of Mode 2",...
        "Root Locus of Mode 3",...
        % "Root Locus of Mode 4",...
        "Imag Axis"
        });
set(leg,"fontsize",8,"location","northwest");
xlb =xlabel("k");
ylb =ylabel("Imp");
set([xlb,ylb],"fontsize",8);

% Compute the flutter speed and which mode corresponding to flutter mode
for imode = 1:neig
    % Find if the real part of phat >0
    loc_list = find( real(pconv(imode,:))>0 );
    % If exist, the number of coloum will no longer be 0
    if size(loc_list,2) > 0
        fprintf("Find the fultter mode is at mode  %.d \n",imode);
        % Corresponding the first loctation is flutter speed
        loc = loc_list(1);
        p_flutter = pconv(imode,loc)
        % loc1 = 0;
            
    end
end

if loc !=0
    uflutter = uvec(loc);
    fprintf("\nThe flutter speed is %.3f m/s \n",uflutter);
    % Retrive the corresponding laplacian non-dimensional frequency
    % The imag part is the  frequency omega        
    f_flutter = imag(p_flutter) * uflutter/ (2*pi*Qip.bref);
    fprintf("\nThe flutter frequency(imag part) is %.3f rad/s \n",f_flutter);
else 
    uflutter = 0;
    fprintf("\nDid not find the flutter speed yet! \n");
    fprintf("\nMaximum Speed now is %.2f m/s\n",uvec(end));
end



[urev,zrev] = reversal(K,Qip,f,CRv,CRd);
fprintf("\nReversal Speed  is %.2f m/s\n",urev);

[udiv,zdiv] = divergence(K,Qip);
fprintf("\nDivergence Speed  is %.2f m/s\n",udiv);




