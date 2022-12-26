function valid_plot(v)
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
    

    t = c * 0.01;
    r = 2.6723;
    Mx = 1.01*E * Ixx .* w2';
    Z = Mx /(r) ;
    fprintf("Maximum Bending moment %.2f \n",max(Mx));
    fprintf("Maximum Shear Force %.2f \n",max(Z));
    sigma = 0.5*t*Mx /Ixx;
   

    figure(50);
    
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperUnits', 'inches');
    x0=15;y0=75;width=600;height=400;
    set(gcf,'units','points','position',[x0,y0,width,height])
    subplot(2,1,1);
    hold on 
    plot(ye,Z,"go-","linewidth",1.8);
    title("Shear Force Distribution","fontsize",15)
    ylabel("Shear Force (N)","fontsize",8)
    a = get(gca,'XTickLabel');
    b = get(gca,'YTickLabel');
    set(gca,'XTickLabel',a,'fontsize',15)
    set(gca,'YTickLabel',b,'fontsize',15)

    subplot(2,1,2);
    plot(ye,Mx,"bs-","linewidth",1.8)
    title("Bending Moment Distribution","fontsize",15)
    ylabel("Moment (Nm)","fontsize",8)
    xlabel('Span coordinate [m]',"fontsize",15)

    a = get(gca,'XTickLabel');
    b = get(gca,'YTickLabel');
    set(gca,'XTickLabel',a,'fontsize',15)
    set(gca,'YTickLabel',b,'fontsize',15)

    print -djpg StressValid.jpg

    Ultimate_norm = 633E6 ;
    factor_normal = Ultimate_norm / max(abs(sigma));
    fprintf("Minimum Safety factor of normal stress is %.2f\n", factor_normal)
end