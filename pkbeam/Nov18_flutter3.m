
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
% If now we compute the uf as divergence 
% The flutter probelm has one solution for divergence
[udiv,zdiv] = divergence(K,Qip)
uf = udiv
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



% To debug programme , use pk_bound plot the k
u = uf;
% how many eigenvalue you used here
neig = 4;
[kbounds] = pk_bounds(u,M,K,Qip,neig)





figure(1);
% PLOT the k VS imag(phat) sorted  
plot(Qip.ktab,imag(pmx),".");
hold on 
plot(kbounds,kbounds,"r*","markersize",8);
plot([0 2],[0 2]);
axis([0 2 0 2]);

%% Explain spy(Qip.Qtab(:,:,1))
figure(2)
spy(Qip.Qtab(:,:,1))

% Just aerppdynamic matrix 
figure(3)
spy(Qip.Qtab(:,:,2))


% dlm : detailed aerodynamic model
load dlmqip
% if you have 24 elements of FEM, 
% you can compare the result with dlmpiq 
% dont forget the constrain for Qdlm.Qtab(75*75)
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

%%%%%% Reduce computation efforts
% Reduced size 
[V,D] = eig(K,M);
ev = diag(D);
nmode = 6;
% New basis

Zm = V(:,1:nmode);
Km = Zm'*K*Zm;
Mm = Zm'* M * Zm;
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
        phat2 = eig(qdyn * Zm'*Q*Zm- Km, Mm.*(uf/b)^2);
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