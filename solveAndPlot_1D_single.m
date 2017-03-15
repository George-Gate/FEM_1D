% Finite Element Method/Finite Difference Method Solver

%% parameters
b=1;
c=0;
k=1;
f=@(x)x.^k;
epsilon=1e-3;
n=2^4;
% Differential Format: central, forward, backward or FEM
dFmtList={'central','forward','backward','FEM'};
dFmt=dFmtList{4};
meshType='shishkin';


%% analytical solution
% depends on epsilon, b, c and k
% get @(x)anaSol(x)
getAnaSol;

%% numerical solution    
% the following depends on dFmt, f(x) and n
% get the coefficient matrices S, C, M and vecf
meshWidth=min(0.49,epsilon/b*2.5*log(n));
getCoeffs;

% the following depends on n, epsilon, b and c
H=epsilon*S+b*C+c*M;

% solve
u=H\vecf;


%% plot
figure();
plot([0;xList;1],[0;u;0],'-o');hold on;

[ax,~,~]=plotyy(0:0.0001:1,real(anaSol(0:0.0001:1)),...
                xList,abs( u-real(anaSol(xList)) ),...
                @(x,y)plot(x,y,'linewidth',2,'color','red'),...
                @(x,y)semilogy(x,y,'g*'));hold off;box on;

% refine plot
legend({'Numerical Solution','Analytical Solution','Error'},'Location','northwest');
title(['$$N=',num2str(N),'\quad \varepsilon=$$',num2str(epsilon),'$$\quad b=',num2str(b),'\quad c=',num2str(c),'\quad f(x)=x^k, k=',num2str(k),'$$  dFmt=',dFmt],'interpreter','latex');
xlabel('$$x$$','interpreter','latex');
ylabel(ax(1),'$$u(x)$$','interpreter','latex','color','black');
ylabel(ax(2),'Error','color','black');
set(ax(1),'fontsize',12,'ylim',[0,1],'Ycolor','black');
set(ax(2),'fontsize',12,'Ycolor','black');