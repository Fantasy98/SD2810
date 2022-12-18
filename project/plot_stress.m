function plot_stress(v,l,b,t,E,G,c)
#function for plotting the normal stress distribution
# Input: 
        #v : deformation vector
        #l : length of beam

    % b = 0.175;
    % t = 0.004; 
    % E = 31.5E9;
    % G = 5.52E9;

    Ixx = 1.4117e+05 * 10^(-12);
    
    ndof = length(v);
    nnode = fix(ndof/3);
    nshape = [1, nnode];
    nelem = nnode-1
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
    

    t = 0.11;
    Mx = E * Ixx .* w2';
    sigma = 0.5*t*Mx /Ixx;
   

    figure(50);
    subplot(2,1,1);
    hold on 
    plot(ye,sigma,"ro-","linewidth",1.5);
    title("Normal Stress Distribution","fontsize",15)
    ylabel("Normal Stress (pa)","fontsize",8)


    subplot(2,1,2);
    plot(ye,Mx,"bs-","linewidth",1.5)
    title("Bending Moment Distribution","fontsize",15)
    ylabel("Moment (N*m)","fontsize",8)

    
    K = 1.45* 0.0212 * c^3 * 6e-4; % torsion constant of beam section [m^4], assume the thickness is t = 6mm
    GK = G*K;
    thetap = diff(theta);
    Mt = GK * thetap;
    Tx = Mt/(0.5*t);
    A = 2*80*110*10^(-6);
    tau = Mt/(2*t*A);

    figure(51);
    subplot(2,1,1);
    hold on 
    plot(ye,tau,"ro-","linewidth",1.5);
    title("Shear Stress Distribution","fontsize",15)
    ylabel("Normal Stress (pa)","fontsize",8)


    subplot(2,1,2);
    plot(ye,Tx,"bs-","linewidth",1.5)
    title("Shear Force Distribution","fontsize",15)
    ylabel("Shear Force (N*m)","fontsize",8)

    % figure(52)
    % plotmode_multi(v);

    % Minimum Safety Factor = Ultimate_stress/ max_Allowed_stress
    % In our case, since the material property is similar to Aluminum, we can assume 
    % Assume the ultimate stress is identical to alumn's 
    Ultimate = 633E6;
    factor_normal = Ultimate / max(abs(sigma));
    factor_shear = Ultimate / max(abs(tau(2:end)));
    fprintf("Minimum Safety factor of normal stress is %.2f\n", factor_normal)
    fprintf("Minimum Safety factor of shear stress is %.2f\n", factor_shear)

end