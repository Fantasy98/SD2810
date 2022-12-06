
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
mhinge = 28/1000;
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

% set deflection of alieron 
% The certia, set 
del = 1:5;
delta = del.*pi/180; % rad
% Set speed
iu = 0; 

for u = 12:.1:20
    iu = iu+1;
    
    q = 0.5* Qip.rho * u * u;
    A = Qip.Qtab(:,:,1);
    for d=del
        % Solve for the deformation
        v = (K-q*A)\(q*f*delta(d));
            

            % Plot deformation 
        v_all(iu,:) = v;
        CR(d,iu) = CRv * v + CRd*d;
        deform_tip(d,iu) = v(end-2);
    end
    uv(iu) = u;
end

% A_eig = A - f*CRv./CRd;
% [V,D] = eig(A_eig , K);
% qrev = 1/max(diag(D));
% urev = sqrt(2*qrev/Qip.rho);


location = abs(deform_tip(1,:)-deform_tip(2,:)) <= 0.001;
urev = uv(location);
% Plot the deformation VS velocity 
figure(1)

for d = del
plot(uv,deform_tip(d,:),"-","linewidth",3.5);
hold on 

end
plot([urev urev],[min(min(deform_tip))-0.2 max(max(deform_tip))+0.2],"r-.","linewidth",3);
grid()
leg = legend({  ["{\delta} ="   num2str(del(1)) "deg"], ...
                ["{\delta} =" num2str(del(2)) "deg"], ...
                ["{\delta} =" num2str(del(3)) "deg"], ...
                ["{\delta} =" num2str(del(4)) "deg"], ...
                ["{\delta} =" num2str(del(5)) "deg"], ...
                ["Reversal Speed"],...
            


            });
axis([min(uv) max(uv) min(min(deform_tip)) max(max(deform_tip))])
set(leg,"location","southwest","fontsize",20,'interpreter', 'tex');
xlab =xlabel("Speed (m/s)","fontsize",18);
ylab =ylabel("Deformation (m)","fontsize",18);
set([xlab ylab],"fontsize",18)
a = get(gca,'XTickLabel');
b = get(gca,'YTickLabel');
set(gca,'XTickLabel',a,'fontsize',20)
set(gca,'YTickLabel',b,'fontsize',20)
print -djpg -r300 Reversal_vs_speed.jpg

A_eig = A - f*CRv./CRd;
[V,D] = eig(A_eig , K);
qrev = 1/max(diag(D));
urev = sqrt(2*qrev/Qip.rho);
fprintf("The reverse velocity is %.2f\n",urev);
% figure(2)
% plot(uv,CR(1,:));
% hold on
% plot([urev urev],[(min(CR(1,:))) max(CR(1,:))],"r-.","linewidth",2);

figure(12)
plotmode_multi(v_all(68,:));
hold on
plotmode_multi(v_all(location,:));
plotmode_multi(v_all(72,:));
grid()
axis([0 l+0.1 -0.1 0.05]);

leg = legend({  ["u = " num2str(uv(68)) "m/s"], ...
                ["reversal =" num2str(uv(70)) "m/s"],....
                ["u = " num2str(uv(72)) "m/s"],....
            
            })
            
set(leg,'location',"southwest","fontsize",18);
a = get(gca,'XTickLabel');
b = get(gca,'YTickLabel');
set(gca,'XTickLabel',a,'fontsize',20)
set(gca,'YTickLabel',b,'fontsize',20)
print -djpg -r300 Reversal_deflection.jpg

return;
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
