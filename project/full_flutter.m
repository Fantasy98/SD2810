clear all

nelem =24;
nnodes = nelem + 1;

l =12.1; % m

c = 0.64;

b = c/2 ; % m

ba = c*0.225/2; % m

mhinge = 5*28e-3;
 
t = 6e-3;

rhop = 1760;

G = 8600E6;
E = 23.9E9;

npmass = 7 ;
m1 = 3*(40.33+6.39+2)/1000; %kg
m2 = 3*(20.06+2*2)/1000; %kg
m3 = 3*(40.33+2*6.39+2*2)/1000; %kg

m_fuselage = 260 - 88;
m_pilot = 53; 
m_parachute = 5;
m_fuse = m_fuselage + m_pilot + m_parachute;
x_coord = (b-ba);
dpm = zeros(npmass,3);
dpm(1,:) = [m1 x_coord 0];dpm(2,:) =[m2 x_coord 0.5*l/6];dpm(3,:) =[m3 x_coord  2*0.5*l/6];
dpm(4,:) =[m2 x_coord  3*0.5*l/6];dpm(5,:) =[m3 x_coord  4*0.5*l/6];dpm(6,:) =[m2 x_coord 5*0.5*l/6];
dpm(7,:) =[m1 x_coord 0.5*l];
dpm(8,:) =[m_fuse 0 l/2];
dpm(9,:)  =[m2 x_coord  (l/12 + l/2)];dpm(10,:) =[m3 x_coord  (2*l/12 + l/2)];dpm(11,:) =[m2 x_coord  (3*l/12 +l/2)];
dpm(12,:) =[m3 x_coord  (4*l/12 + l/2)];dpm(13,:) =[m2 x_coord  (5*l/12 + l/2)];dpm(14,:) =[m1 x_coord  l];

ndof = 3*nnodes;
B = []

% retrieve system matrices
[M,K,Z,Qip,f,CRv,CRd] = nwing(B, l, b, t, ba, mhinge, rhop, E, G, nelem, dpm);


e1 = zeros(ndof,1);e1(1:3:end) = 1;
e3 = zeros(ndof,1);e3(3:3:end) = 1;

## Assemble constrain matrix
yn = 0.5*l*linspace(-1,1,nnodes);idof = 3*(nelem/2)+3;B = zeros(3,ndof);
B(1,:) = e1';
B(2,2:3:end) = 1;
B(2,1:3:end) = yn;
B(3,:) = e3';

Z = null(B);
nz=7;u= 350/3.6;

k=0;
Stail= 0.86;le = 2.15;
bt = 0.5 * Stail/le;

for i=1:length(Qip.ktab)
    kw = Qip.ktab(i);
    Qtail = tail_matrix(kw,nelem);
    Q = (Qip.Qtab(:,:,i) + Qtail);
    Qall.Qtab(:,:,i) = Z'*Q*Z;
    ke = kw*0.32/bt;
    Qall.ktab(i) = kw;
end
Qall.rho = Qip.rho;
Qall.bref = Qip.bref;


nmode = 12;
[Km,Mm,Zm,mQip]=ReduceDim(Z'*M*Z,Z'*K*Z,Qall,nmode);
neig= 6
[ucrit,pcrit,zcrit,pconv,uvec] = flutter(Mm,Km,mQip,neig,iter = 100)
% [Q,vtot,ve,alfa0,delta0] = flight_load(K,M,B,Qip,f,nelem,u,nz,k);
uvec = uvec *3.6;
figure(1)
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
x0=15;y0=75;width=600;height=400;
set(gcf,'units','points','position',[x0,y0,width,height])

plot(uvec,real(pconv),"-o","markersize",7.5)
hold on 
plot([min(uvec) max(uvec)+40],[0 0],"k-.","linewidth",1.5)
plot([350 350],[-0.3,0.05],"r-.","linewidth",2)
% plot([70/3.6 70/3.6],[-0.3,0.05],"r-.","linewidth",3)
axis([min(uvec) max(uvec)+30 -0.06 0.01])
leg = legend({
                "Mode 1",
                "Mode 2",
                "Mode 3",
                "Mode 4",
                "Mode 5",
                "Mode 6",
                "Real=0",
                "Diving Speed"
})
set(leg,"fontsize",3,"location","southeast")
xlab = xlabel("Speed (km/h)")
ylab = ylabel("Real p")
set([xlab,ylab],"fontsize",12);
a = get(gca,'XTickLabel');
b = get(gca,'YTickLabel');
set(gca,'XTickLabel',a,'fontsize',15)
set(gca,'YTickLabel',b,'fontsize',15)

print -djpg flutterfree.jpg

figure(2)
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
x0=15;y0=75;width=600;height=400;
set(gcf,'units','points','position',[x0,y0,width,height])

plot(real(pconv),imag(pconv),"o","markersize",7.5)
hold on 
plot([0 0],[-0.1 0.1],"k-.","linewidth",1.5)
% plot([350 350],[-0.3,0.05],"r-.","linewidth",2)
% plot([70/3.6 70/3.6],[-0.3,0.05],"r-.","linewidth",3)
% axis([min(uvec) max(uvec)+30 -0.06 0.01])
leg = legend({
                "Mode 1",
                "Mode 2",
                "Mode 3",
                "Mode 4",
                "Mode 5",
                "Mode 6",
                "Imag Axis",
                
})
set(leg,"fontsize",3,"location","southeast")
xlab = xlabel("Imag p")
ylab = ylabel("Real p")
set([xlab,ylab],"fontsize",12);
a = get(gca,'XTickLabel');
b = get(gca,'YTickLabel');
set(gca,'XTickLabel',a,'fontsize',15)
set(gca,'YTickLabel',b,'fontsize',15)

print -djpg flutterfree2.jpg