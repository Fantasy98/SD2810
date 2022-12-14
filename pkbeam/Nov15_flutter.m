
% Nov 15: Flutter 
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
% #########################################
% Flutter case
% Plot how do eigenvalue change when k= 0 : 2 since max k in Qip.ktip is 2
% set speed



uf = 15;
for ik = 1:length(Qip.ktab)
    k= Qip.ktab(ik);
    qdyn = 0.5*Qip.rho*uf *uf;
    phat2 = eig(qdyn*Qip.Qtab(:,:,ik)-K, M.*(uf/b)^2);
    phat = sqrt(phat2);
    phat = phat .* sign(imag(phat));

    % Sort the phat by ascending order and store the order of index
    [psort ipsort] = sort(imag(phat));
    % Store it into the pmx
    pmx(:,ik) = phat(ipsort);
end

figure(1);
% PLOT the k VS imag(phat) sorted  
plot(Qip.ktab,imag(pmx),".");
hold on 
plot([0 2],[0 2]);
axis([0 2 0 2])
return;
ev = eig(K,M);
omega = sqrt(ev)/(2*pi);
uf = 15;
for iu = 1:20

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
        
        uvec(iu) = uf;
        
        end 
    
    
    end
    uf = uf+0.5;
end

% Plot the root loot of each phat
figure(2);
for imode = 1:3
plot( 
    % [0,real(pconv(imode,:))],...
      [omega(imode).*b/uvec(1), imag(pconv(imode,:))],"o-","linewidth",1.5);

hold on 
end
plot([0 0],[0 1],"k-","linewidth",1.5);
axis([-0.5 0.5 0 0.5]);
leg = legend({
        "Root Loot of Mode 1",...
        "Root Loot of Mode 2",...
        "Root Loot of Mode 3",...
        "Imag Axis"
        });
set(leg,"fontsize",16);
xlabel("k");
ylabel("Imp");








