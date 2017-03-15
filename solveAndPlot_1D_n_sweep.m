% Finite Element Method/Finite Difference Method Solver
% Differential Format: central, forward, backward or FEM
% dFmtList={'central','forward','backward','FEM'};
dFmt='FEM';
meshType='shishkin';
extraPointPos=0.9;
%% parameters
b=1;
c=0;
k=2;
f=@(x)x.^k;
epsilon=1e-6;

%% analytical solution
% depends on epsilon, b, c and k
% get @(x)anaSol(x)
getAnaSol;

%% n - sweep
nList=2.^(5)';
Err=zeros(size(nList));
numSol=cell(size(nList));
legendList=cell(length(nList)+1,1);

for i=1:length(nList)
    n=nList(i);
    meshWidth=min(0.49,epsilon/b*2.5*log(n));
    % the following depends on dFmt and n
    % get the coefficient matrices S, C, M and vecf
    getCoeffs;
    
    % the following depends on n, epsilon, b and c
    H=epsilon*S+b*C+c*M;


    % solve
    % tic;
    u=H\vecf;
    % toc;
    
    % record
    numSol{i}.n=n;
    numSol{i}.u=u;
    numSol{i}.xList=xList;
    legendList{i}=['N=',num2str(N)];
    
end

%% plot
figure('position',[100 100 940 360]);
for i=1:length(nList)
    plot(numSol{i}.xList,anaSol(numSol{i}.xList)-numSol{i}.u,'-o');
end
xlabel('$$x$$','interpreter','latex');ylabel('anaSol - numSol');
title('Error between numerical and analytical solution.');

figure('position',[0 0 940 360]);
for i=1:length(nList)
    plot([0;numSol{i}.xList;1],[0;numSol{i}.u;0],'-o');hold on;
end
plot(0:1e-6:1,anaSol(0:1e-6:1),'linewidth',2);
legendList{end}='Analytical Solution';
if strcmp(dFmt,'FEM')
    xlabel('$$x$$','interpreter','latex');ylabel('$$u^{(\mathrm{FEM})}$$','interpreter','latex');
else
    xlabel('$$x$$','interpreter','latex');ylabel('$$u^{(\mathrm{FDM})}$$','interpreter','latex');
end
title(['$$\varepsilon=\mathrm{',num2str(epsilon,'%1.1E'),'}\quad b=',num2str(b),'\quad c=',num2str(c),'\quad k=',num2str(k),'$$  \quad dFmt=',dFmt,' \quad meshType=',meshType],'interpreter','latex');
set(gca,'fontsize',12);
legend(legendList,'location','north');

