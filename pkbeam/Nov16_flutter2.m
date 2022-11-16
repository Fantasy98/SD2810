
% Nov 15: Flutter 
clear all;

% setup geometry and structural properties
% number of finite elements requested should be a multiple of 3
nelem = 10;
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
% #########################################
% Flutter case
% Plot how do eigenvalue change when k= 0 : 2 since max k in Qip.ktip is 2
% set speed



% uf = 15;
% for ik = 1:length(Qip.ktab)
%     k= Qip.ktab(ik);
%     qdyn = 0.5*Qip.rho*uf *uf;
%     phat2 = eig(qdyn*Qip.Qtab(:,:,ik)-K, M.*(uf/b)^2);
%     phat = sqrt(phat2);
%     phat = phat .* sign(imag(phat));

%     % Sort the phat by ascending order and store the order of index
%     [psort ipsort] = sort(imag(phat));
%     % Store it into the pmx
%     pmx(:,ik) = phat(ipsort);
% end

% figure(1);
% % PLOT the k VS imag(phat) sorted  
% plot(Qip.ktab,imag(pmx),".");
% hold on 
% plot([0 2],[0 2]);
% axis([0 2 0 2])

ev = eig(K,M);
omega = sqrt(ev)/(2*pi);
uf = 15;
for iu = 1:100

    % At each speed, compute the dynamic pressure
    for imode = 1:4
    % Compute dynamic pressure
    qdyn = 0.5*Qip.rho*uf *uf;
    % Reduced freq definition: k = omega * b/u
    k = omega(imode)*b/uf;
    
        for iter = 1:10
        % Interploation by given reduced frequency to derive aerodynamic force Q
        Q = ipolQk(Qip,k);
        % Solve eigenvalue problem to get non-dimensional laplace coefficient
        phat2 = eig(qdyn*Q-K, M.*(uf/b)^2);
        % the squre of phat is complex number
        % Real part indicates the stability
        % imag part indicates the oscillation frequency
        phat = sqrt(phat2);
        % pass the sign of imag part then it's easier to sort them
        phat = phat .* sign(imag(phat));

        % Sort the phat by ascending order 
        % and store the order of index
        [psort ipsort] = sort(imag(phat));
        % set k as the lower boundary for this iteration
        k= imag(phat(ipsort(imode)));
        % Store the phat
        pconv(imode,iu) = phat(ipsort(imode));
            
        judge = real(pconv(:,iu)) < 0;
            if length(unique(judge)) > 1
                uflutter = uf;
            end
        
        
        
        end 
    
    
    end
    uvec(iu) = uf;
    uf = uf+0.1;
    
end

% Plot the root loot of each phat
figure(2);
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


% Plot the real part of phat versus velocity
figure(3)
for imode = 1:4
    plot(uvec, real(pconv(imode,:)),"o-","linewidth",0.8,"markersize",4.5);
    hold on 
    end
plot([14 25],[0,0],"k-.","linewidth",1.5);
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



