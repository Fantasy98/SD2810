function plot_stress(v)
#function for plotting the normal stress distribution
# Input: 
        #v : deformation vector
        #l : length of beam

    % b = 0.175;
    % t = 0.004; 
    % E = 31.5E9;
    % G = 5.52E9;
    l =12.1; % m

    c = 0.64;
# According to the technic report, the main wing is assumed as ellipse with chord = 500mm
    b = c/2 ; % m
# The alieron span is defined as c*E, E = 0.225 in W12C-JARA-0001.pdf P
    t = 6e-3;
    G = 8600E6;
    E = 23.9E9;
    
    Ixx = 1.6117e+05 * 10^(-12);
    Ixx = 1.6E-5;
    
    ndof = length(v);
    nnode = fix(ndof/3);
    nshape = [1, nnode];
    nelem = nnode-1;
    # We need elemental deformation! See the equation (59) in the LectureNote


   
    % Extract bending and torsion dofs section by section
    w = zeros(nshape);
    theta = zeros(nshape);
    for k = 1:nnode
        w(k) = v(1+3*(k-1)); % nodal deflection
        w1(k) = v(2+3*(k-1)); % nodal deflection curvature

        
        theta(k) = 180/pi * v(3+3*(k-1)); % nodal twist
        % theta(k) = v(3+3*(k-1)); % nodal twist
    end
    % Compute normal stress along span
    % Moment = -E * I * w'',
    % w'' is the curvature of beam
    % Then the normal stress can be derived as sigma = M * z/I
    % z = 0.5 * t if we only consider the outer part
  
    
    yp = linspace(0.0, l, nnode);
    ye = linspace(0,l,nelem);

    # Interpolation of w'' 
    le = yp(2)-yp(1);
    w2 = diff(w1)/le; 
    t = c * 0.17;
  
    Mx = E * Ixx .* w2';
  
    sigma = 0.5*t.*Mx /Ixx;
   

    figure(50);
    
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperUnits', 'inches');
    x0=15;y0=75;width=600;height=400;
    set(gcf,'units','points','position',[x0,y0,width,height])
    subplot(2,1,1);
    hold on 
    plot(ye,sigma,"ro","linewidth",1.8,"markersize",6,'MarkerFaceColor','r');
    title("Normal Stress Distribution","fontsize",15)
    ylabel("Normal Stress [pa]","fontsize",8)
    a = get(gca,'XTickLabel');
    b = get(gca,'YTickLabel');
    set(gca,'XTickLabel',a,'fontsize',15)
    set(gca,'YTickLabel',b,'fontsize',15)

    
    K = 1.45* 0.0212 * c^3 * 0.6e-3; % torsion constant of beam section [m^4], assume the thickness is t = 6mm
    GK =G*K;


    thetap = (1/le).*diff(theta)*pi/180;
    Mt = GK * thetap ;
    
    % A = 54000*0.55*0.3*10^(-9);
    A = 54000*t*0.6*10^(-9);
    tau = Mt/(2*A);
    

    subplot(2,1,2);
    plot(ye,tau,"bs","linewidth",1.8,"markersize",6,'MarkerFaceColor','blue')
    title("Shear Stress Distribution","fontsize",15)
    ylabel("Shear Stress [pa]","fontsize",8)
    xlabel('Span coordinate [m]',"fontsize",15)
    a = get(gca,'XTickLabel');
    b = get(gca,'YTickLabel');
    set(gca,'XTickLabel',a,'fontsize',15)
    set(gca,'YTickLabel',b,'fontsize',15)

    % print -djpg Stress73.jpg
    % print -djpg Stress53.jpg


    % Minimum Safety Factor = Ultimate_stress/ max_Allowed_stress
    % In our case, since the material property is similar to Aluminum, we can assume 
    % Assume the ultimate stress is identical to alumn's 

    % Ultimate_norm = 633E6
    % Ultimate_norm = 53007 * 0.5*t/Ixx
    Ultimate_norm = 208E6
   
    Ultimate_shear = 158E6
    % Ultimate_shear = 31E6
    
    factor_normal = Ultimate_norm / max(abs(sigma));
    factor_shear = Ultimate_shear / max(abs(tau));
    fprintf("Minimum Safety factor of normal stress is %.2f\n", factor_normal)
    fprintf("Minimum Safety factor of shear stress is %.2f\n", factor_shear)

end