function plot_stress(v,l)
#function for plotting the normal stress distribution
# Input: 
        #v : deformation vector
        #l : length of beam

    b = 0.175;
    t = 0.004; 
    E = 31.5E9;
    G = 5.52E9;


    
    ndof = length(v);
    nnode = fix(ndof/3);
    nshape = [1, nnode];

    % Extract bending and torsion dofs section by section
    w = zeros(nshape);
    theta = zeros(nshape);
    for k = 1:nnode
        w(k) = v(1+3*(k-1)); % nodal deflection
        w2(k) = v(2+3*(k-1)); % nodal deflection curvature
        
        theta(k) = 180/pi * v(3+3*(k-1)); % nodal twist
    end

    % Compute normal stress along span
    % Moment = -E * I * w'',
    % w'' is the curvature of beam
    % Then the normal stress can be derived as sigma = M * z/I
    % z = 0.5 * t if we only consider the outer part
    yp = linspace(0.0, l, nnode);
    Ixx= (2*b*t^3)/12;


    % delta_theory = P(end-2)*l^3 /(3*E*Ixx)

    Mx = -E * Ixx .* w2 .*yp;
    sigma = Mx .*0.5*t./Ixx;
    
    % Compute shear stress along span
    % theta = T * L / (J*G)
    beta = 1/3;
    alpha = 1/3;
    Tx = theta*beta * 2*b * t^3 * G ./yp;
    tau = Tx./(alpha * (2*b) * t^2);
    clf;
    figure(50);
    subplot(2,1,1);
    hold on 
    plot(yp,sigma(end:-1:1),"ro-","linewidth",1.5);
    title("Normal Stress Distribution","fontsize",15)
    ylabel("Normal Stress (pa)","fontsize",8)


    subplot(2,1,2);
    plot(yp,Mx(end:-1:1),"bs-","linewidth",1.5)
    title("Bending Moment Distribution","fontsize",15)
    ylabel("Moment (N*m)","fontsize",8)


    figure(51);
    subplot(2,1,1);
    hold on 
    plot(yp,tau(end:-1:1),"ro-","linewidth",1.5);
    title("Shear Stress Distribution","fontsize",15)
    ylabel("Normal Stress (pa)","fontsize",8)


    subplot(2,1,2);
    plot(yp,Tx(end:-1:1),"bs-","linewidth",1.5)
    title("Shear Force Distribution","fontsize",15)
    ylabel("Shear Force (N*m)","fontsize",8)

    figure(52)
    plotmode(v);

    % Minimum Safety Factor = Ultimate_stress/ max_Allowed_stress
    % In our case, since the material property is similar to Aluminum, we can assume 
    % Assume the ultimate stress is identical to alumn's 
    Ultimate = 90e6;
    factor_normal = Ultimate / max(abs(sigma));
    factor_shear = Ultimate / max(abs(tau(2:end)));
    fprintf("Minimum Safety factor of normal stress is %.2f\n", factor_normal)
    fprintf("Minimum Safety factor of shear stress is %.2f\n", factor_shear)

end