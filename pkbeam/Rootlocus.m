function Rootlocus(pconv,neig)
% Function for plot root locus
% Input: 
%   pconv: laplacian frequency computed by pk_bisection 
%   neig:  Number of modes to be ploted
% Output:
%   Figure of root locus
    x0=100;y0=850;width=550;height=400;
    set(gcf,'position',[x0,y0,width,height]);
    for i = 1:neig
    plt(i)=plot(real(pconv(i,:)),imag(pconv(i,:)),"o-","linewidth",1,"markersize",7.5,"Displayname",["Mode" num2str(i)]);

    hold on 
    end
    plt(neig+1)=plot([0 0],[0 1],"k-","linewidth",3.5,"Displayname","Imag Axis");
    axis([-0.2 0.2 0 0.5]);
    leg = legend;

    set(leg,"fontsize",20,"location","northwest");
    xlb =xlabel("Rep");
    ylb =ylabel("Imp");
    set([xlb,ylb],"fontsize",18);
    a = get(gca,'XTickLabel');
    b = get(gca,'YTickLabel');
    set(gca,'XTickLabel',a,'fontsize',20)
    set(gca,'YTickLabel',b,'fontsize',20)
    
end