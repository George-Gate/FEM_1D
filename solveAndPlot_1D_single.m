% Finite Element Method/Finite Difference Method Solver

%% parameters
b=0;
c=1;
k=0;
f=@(x)x.^k;
epsilon=1e-5;
n=2^8;
% Differential Format: central, forward, backward or FEM
dFmtList={'central','forward','backward','FEM'};
dFmt=dFmtList{4};
meshType='uniform';


%% analytical solution
% depends on epsilon, b, c and k
% get @(x)anaSol(x)
getAnaSol;

%% numerical solution    
% the following depends on dFmt, f(x) and n
% get the coefficient matrices S, C, M and vecf
if (b)
    meshWidth=min(0.49,epsilon/b*2.5*log(n));
else
    meshWidth=min(1/3.1,sqrt(epsilon/c)*2.5*log(n));
end
getCoeffs;

% the following depends on n, epsilon, b and c
H=epsilon*S+b*C+c*M;

% solve
u=H\vecf;


%% plot
figure();
plot([0;xList;1],[0;u;0],'-o');hold on;

switch meshType
    case 'shishkin'
        w=-epsilon*log(epsilon);
        xSample=[linspace(0,1-w,5*nList(end))';linspace(1-w,1,5*nList(end))'];
    case '2sideShishkin'
        w=-sqrt(epsilon)*log(epsilon)/2;
        xSample=[linspace(0,w,5*nList(end))';linspace(w,1-w,5*nList(end))';linspace(1-w,1,5*nList(end))'];
    otherwise
        xSample=linspace(0,1,5*nList(end))';
end
[ax,~,~]=plotyy(xSample,real(anaSol(xSample)),...
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