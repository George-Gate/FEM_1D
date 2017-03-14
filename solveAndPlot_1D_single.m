% Finite Element Method/Finite Difference Method Solver

%% parameters
b=0;
c=1;
k=0;
f=@(x)x.^k;
epsilon=1e-4;
n=2^6;
% Differential Format: central, forward, backward or FEM
dFmtList={'central','forward','backward','FEM'};
dFmt=dFmtList{3};


%% analytical solution
% depends on epsilon, b, c and k
% get @(x)anaSol(x)
getAnaSol;

%% numerical solution    
% the following depends on dFmt, f(x) and n
h=1/n;
% get the coefficient matrices S, C, M and vecf
getCoeffs;

% the following depends on n, epsilon, b and c
H=epsilon*S+b*C+c*M;

% solve
u=H\vecf;


%% plot
figure();
plot(0:h:1,[0;u;0],'-o');hold on;

[ax,~,~]=plotyy(0:0.0001:1,real(anaSol(0:0.0001:1)),...
                h:h:1-h,abs( u'-real(anaSol(h:h:1-h)) ),...
                @(x,y)plot(x,y,'linewidth',2,'color','red'),...
                @(x,y)semilogy(x,y,'g*'));hold off;box on;

% refine plot
legend({'Numerical Solution','Analytical Solution','Error'},'Location','northwest');
title(['$$n=',num2str(n),'\quad \varepsilon=$$',num2str(epsilon),'$$\quad b=',num2str(b),'\quad c=',num2str(c),'\quad f(x)=x^k, k=',num2str(k),'$$  dFmt=',dFmt],'interpreter','latex');
xlabel('$$x$$','interpreter','latex');
ylabel(ax(1),'$$u(x)$$','interpreter','latex','color','black');
ylabel(ax(2),'Error','color','black');
set(ax(1),'fontsize',12,'ylim',[0,1],'Ycolor','black');
set(ax(2),'fontsize',12,'Ycolor','black');