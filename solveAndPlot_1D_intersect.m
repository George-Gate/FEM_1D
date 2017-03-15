% Finite Element Method/Finite Difference Method Solver
% Differential Format: central, forward, backward or FEM
dFmtList={'central','forward','backward','FEM'};
meshTypeList={'uniform','uniformP1','uniformP2','shishkin','2sideShishkin'};
dFmt=dFmtList{4};

extraPointPos=0.9;
%% parameters
b=0;
c=1;
k=2;
f=@(x)x.^k;
epsilon=1e-12;
n=2^5;

%% analytical solution
% depends on epsilon, b, c and k
% get @(x)anaSol(x)
getAnaSol;

%% solve
numSol=cell(2,1);topSol=1;
legendList=cell(length(numSol)+1,1);

for i=[1,2]
    meshType=meshTypeList{i};
    %meshWidth=min(0.49,epsilon/b*2.5*log(n));
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
    numSol{topSol}.meshType=meshType;
    numSol{topSol}.u=u;
    numSol{topSol}.xList=xList;
    legendList{topSol}=meshType;
    topSol=topSol+1;
    
end

%% plot
% interporlate
for i=1:topSol-1
    numSol{i}.pp=griddedInterpolant([0;numSol{i}.xList;1],[0;numSol{i}.u;0]);
end
[intX,intY]=linkIntersection([0;numSol{2}.xList;1],numSol{1}.pp([0;numSol{2}.xList;1]),[0;numSol{2}.u;0]);
lineInt=griddedInterpolant(intX,intY);
% plot error
figure('position',[0 0 940 700]);
subplot(2,1,1);
xList=numSol{1}.xList(2):0.01/n:numSol{1}.xList(end-1);
for i=1:topSol-1
    plot(xList,anaSol(xList)-numSol{i}.pp(xList));hold on;
end
factor=1e2;  % amplify error of Intersection Line
plot(xList,factor*(anaSol(xList)-lineInt(xList)),'--');
plot(intX(3:end-3),factor*(anaSol(intX(3:end-3))-intY(3:end-3)),'*');
legendList{topSol}=['Intersection x ',num2str(factor,'%.1E')];
legendList{topSol+1}='Intersection Line data point';
legend(legendList,'location','southwest');
xlabel('$$x$$','interpreter','latex');ylabel('anaSol - numSol');
title('Error between numerical and analytical solution.');

% plot solution
subplot(2,1,2);
for i=1:topSol-1
    plot([0;numSol{i}.xList;1],[0;numSol{i}.u;0],'-o');hold on;
end
plot(0:1e-6:1,anaSol(0:1e-6:1),'linewidth',2);
plot(intX,intY,'--*');
legendList{topSol}='Analytical Solution';
legendList{topSol+1}='Intersection Line';
if strcmp(dFmt,'FEM')
    xlabel('$$x$$','interpreter','latex');ylabel('$$u^{(\mathrm{FEM})}$$','interpreter','latex');
else
    xlabel('$$x$$','interpreter','latex');ylabel('$$u^{(\mathrm{FDM})}$$','interpreter','latex');
end
title(['$$\varepsilon=\mathrm{',num2str(epsilon,'%1.1E'),'}\quad b=',num2str(b),'\quad c=',num2str(c),'\quad k=',num2str(k),'$$  \quad dFmt=',dFmt,' \quad n=',num2str(n)],'interpreter','latex');
set(gca,'fontsize',12);
legend(legendList,'location','north');

