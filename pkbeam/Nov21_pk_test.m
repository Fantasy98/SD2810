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
npmass = 7 ;
m1 = (40.33+6.39+2)/1000; %kg
m2 = (20.06+2*2)/1000; %kg
m3 = (40.33+2*6.39+2*2)/1000; %kg
x_coord = (280-175)/1000;
dpm = zeros(npmass,3);
dpm(1,:) = [m1 x_coord 0];
dpm(2,:) =[m2 x_coord 27/100];
dpm(3,:) =[m3 x_coord  53/100];
dpm(4,:) =[m2 x_coord  80/100];
dpm(5,:) =[m2 x_coord  106/100];
dpm(6,:) =[m2 x_coord 133/100];
dpm(7,:) =[m2 x_coord 160/100];

% Test before/after adding hingmass 
% In theory, since hinges are after the elastic axis, the flutter speed should be decreased
dpm = zeros(npmass,3);
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
setenv ("OCTAVE_LATEX_DEBUG_FLAG", "1");
figure(1)
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
% set(gcf, 'PaperPosition', [6 2 6 2]);
plot(Qip.ktab,imag(pmx(1:4,:)),"o-","linewidth",1.2,"markersize",7.5);
hold on 
% Plot the kbounds
plot(kbounds(1:5),kbounds(1:5),"rp","markersize",15,'MarkerFaceColor','red');
% Plot the Imp = k line as reference
plot([0 1.5],[0 1.5],"k","linewidth",1.2);
% Only show k = 0~2 since we only interested in the first 5 modes
leg = legend({["Im p1"]...,
              ["Im p2"]...,
              ["Im p3"]...,
              ["Im p4"]...,
            %   ["Im {{p_5}}"]...,
              ["K boundaries"]...,
              ["Im p(k) =k"]...,

            })
set(leg,"location","northwest","fontsize",18,'Interpreter','Latex');
xlab = xlabel("k","fontsize",15);
ylab = ylabel('Im p',"fontsize",15);
set([xlab ylab],'Interpreter','tex',"fontsize",20)
axis([0 1.5 0 1.5]);
a = get(gca,'XTickLabel');
b = get(gca,'YTickLabel');
set(gca,'XTickLabel',a,'fontsize',20)
set(gca,'YTickLabel',b,'fontsize',20)
print -djpg -r300 Imp_kbound.jpg

return;
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