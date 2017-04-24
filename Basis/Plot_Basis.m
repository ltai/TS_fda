
% just temp to plot basis functions
%
% created on Oct. 14, 2015

close all;
fig = figure;



t = linspace(0,1,101)';

%markers = {'+','o','*','.','x','s','d','^','v','>','<','p','h'};
markers = {'+','o','*','x','s','d','v','>','<','p','h'};
%colors = {'y','m','c','r','g','b','k'};
colors = {'r','g','b','c','m','y','k'};

%  for legender polynomial
My_Figure(8,6);
if(false)
    for i=1:5
        plot(t,Legendre_Poly(:,i),markers{i},'color',colors{i},'linewidth',2); hold on;
    end 

    title('\fontsize{36}Shifted Legendre Basis');
    xlabel('\fontsize{36}t');
    ylabel('\fontsize{36}values');
    legend('\fontsize{22}basis 1','\fontsize{22}basis 2','\fontsize{22}basis 3'...
        ,'\fontsize{22}basis 4','\fontsize{22}basis 5','location','Best'); %'NorthEastOutside'

    ylim([-2,1.5]);
    print(fig,'-dpng','sLegendre_basis');
end


% for Orthogonal exponential basis
My_Figure(8,6); 
for i=1:5
    plot(t,Ortho_Exp_Poly(:,i),markers{i},'color',colors{i},'linewidth',2); hold on;
end 

title('\fontsize{36}Gram-Schmidt to \{e^{a*t}\}');
xlabel('\fontsize{36}t');
ylabel('\fontsize{36}values');
legend('\fontsize{22}basis 1','\fontsize{22}basis 2','\fontsize{22}basis 3'...
    ,'\fontsize{22}basis 4','\fontsize{22}basis 5','location','Best'); %'NorthEastOutside'
ylim([-2.5,6.5]);

print(fig,'-dpng','Exp_basis');
