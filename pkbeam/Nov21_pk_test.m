% Nov 21 testing pk method functions
% pk_bounds and pk_bisect

clear all;

% setup geometry and structural properties
% number of finite elements requested should be a multiple of 3
nelem = 15;
nnodes = nelem + 1;
l =1.6; % m
b = 0.175; % m
ba = 0.03 ;

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

%###############################
%Flutter
% 1 Check pk_bounds function
# Try the first 5 mode using pk_bounds when speed = 15m/s
% First we derive the upper and lower bounds for first 5 modes.
neig = 5;
u = 15;
[kbounds] = pk_bounds(u,M,K,Qip,neig)
Mr = M.*(u/Qip.bref)^2;
% Then We derive all Imp(k) according to existing reduced frequency
for ik = 1: length(Qip.ktab)
    k = Qip.ktab(ik);
    qdyn = 0.5 * Qip.rho* u*u;
    Kr = qdyn*Qip.Qtab(:,:,ik)-K;
    phat2 = eig(Kr,Mr);
    phat = sqrt(phat2);
    phat = phat.*sign(imag(phat));
    % Sort all phat by imag part
    [psort,ipsort] = sort(imag(phat));
    % Store stored result into pmx for each k
    %  
    pmx(:,ik) = phat(ipsort);
end 

% figure(1)
% plot(Qip.ktab,imag(pmx),".");
% hold on 
% % Plot the kbounds
% plot(kbounds,kbounds,"r*","markersize",8.5);
% % Plot the Imp = k line as reference
% plot([0 2],[0 2]);
% % Only show k = 0~2 since we only interested in the first 5 modes
% axis([0 2 0 2]);


%#######################################
% 2 Testing pk_bisect() 
% Note that ieig is the index of eigvalue, 
% which means it can only solve for single mode 

% Remain Question: How to derive the vector? 
u = 13; 
i = 0.5
for iu = 1:15
    
    neig = 4;
    [kbounds] = pk_bounds(u,M,K,Qip,neig);
    for ieig = 1:neig
        [pu,vu] = pk_bisect(ieig,u,M,K,Qip,kbounds);
        pconv(ieig,iu) = pu;
    end

    u = u + i;
    uvec(iu) = u;
end



figure(10);
for imode = 1:4
plot( real(pconv(imode,:)),...
      imag(pconv(imode,:)),"o-","linewidth",0.8,"markersize",4.5);

hold on 
end
plot([0 0],[0 1],"k-","linewidth",1.5);
axis([-0.2 0.2 0 1]);
leg = legend({
        "Root Locus of Mode 1",...
        "Root Locus of Mode 2",...
        "Root Locus of Mode 3",...
        "Root Locus of Mode 4",...
        "Imag Axis"
        });
set(leg,"fontsize",8,"location","northwest");
xlb =xlabel("k");
ylb =ylabel("Imp");
set([xlb,ylb],"fontsize",8);



figure(11);
for imode = 1:4
    plot(uvec, real(pconv(imode,:)),"o-","linewidth",0.8,"markersize",4.5);
    hold on 
    end
plot([10 25],[0,0],"k-.","linewidth",1.5);
leg = legend({
        "Real part of Mode 1",...
        "Real part of Mode 2",...
        "Real part of Mode 3",...
        "Real part of Mode 4",...
        "Re(p) = 0"
        });
set(leg,"fontsize",8,"location","southeast");
xlabel("u (m/s)");
ylabel("Real(p)");
axis([14 25]);