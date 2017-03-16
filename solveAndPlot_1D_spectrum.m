% Finite Element Method/Finite Difference Method Solver

%% parameters
b=0;
c=1;
k=0;
f=@(x)x.^k;
epsilon=1e-12;
n=1;
dFmt='FEM+Spectrum';
meshType='2sideShishkin';
cutOff=[20;5;20;10;5;5];
sigma=1.04;  % width factor of shishkin mesh

%% analytical solution
% depends on epsilon, b, c and k
% get @(x)anaSol(x)
getAnaSol;

%% numerical solution    
% the following depends on dFmt, f(x) and n
% get the coefficient matrices S, C, M and vecf
if b==0
    meshWidth=min(0.49,sigma*sqrt(epsilon)*max(cutOff(1:2*n)));
else
    meshWidth=min(0.49,sigma*epsilon*max(cutOff(1:2*n)));
end
getCoeffs;

% the following depends on n, epsilon, b and c
H=epsilon*S+b*C+c*M;

% solve
tic;
u=H\vecf;
disp(['Time used to solve linear system: ',num2str(toc),'s']);
%% get numSol
tic;
% set sampling points
xList_f=[0;xList;1];
sRate=500;    % how many sampling points should be used in each mesh grid 
gridID=ones(N+1,1); % index of xList(i-1) in xSample
xSample=0;
for i=1:N
    tmp=linspace(xList_f(i),xList_f(i+1),sRate)';
    xSample=[xSample;tmp(2:end)];
    gridID(i+1)=length(xSample);
end
clear tmp;
numSol=zeros(size(xSample));
% linear basis
for i=1:N-1
    i1=gridID(i);i2=gridID(i+1);i3=gridID(i+2);
    numSol(i1  :i2)=numSol(i1  :i2)+u(i)*(xSample(i1  :i2)-xList_f(i))/hList(i);
    numSol(i2+1:i3)=numSol(i2+1:i3)+u(i)*(xList_f(i+2)-xSample(i2+1:i3))/hList(i+1);
end

% Lobatto basis
% calc Lagendre Polynomial first
tmpN=0:max(cutOff(1:N))+1;
if ~exist('legendreMatrix','var') || size(legendreMatrix,2)<length(tmpN) || size(legendreMatrix,1)~=sRate
    % reuse previous result if possible
    tmpX=linspace(-1,1,sRate);
    [tmpN,tmpX]=meshgrid(tmpN,tmpX);
    legendreMatrix=legendreP(tmpN,tmpX);
end
% combine then into numSol
for m=1:N
    i1=gridID(m);i2=gridID(m+1);
    tmpNlist=(1:cutOff(m))';
    
    numSol(i1:i2)=numSol(i1:i2)+sum( (legendreMatrix(:,tmpNlist+2) - legendreMatrix(:,tmpNlist)) * diag(  u(fun2id.psi{m})./sqrt(4*tmpNlist+2)  ),2);
end
clear i1 i2 i3 sRate order tmpX tmpN tmpNlist;

absErrNA=abs( numSol(2:end-1)-real(anaSol(xSample(2:end-1))) );

disp(['Time used to combine the final result: ',num2str(toc),'s']);

%% plot
figure();
plot(xList_f,[0;u(1:N-1);0],'o');hold on;
[ax,~,~]=plotyy(xSample,real(anaSol(xSample)),...
                xSample(2:end-1),absErrNA,...
                @(x,y)plot(x,y,'linewidth',2,'color','red'),...
                @(x,y)semilogy(x,y,'g*'));
plot(xSample,numSol,'k--');hold off;box on;

% refine plot
legend({'Numerical Solution','Analytical Solution','Error'},'Location','northwest');
title({['\centerline{$$N=',num2str(N),'\quad \varepsilon=$$',num2str(epsilon),'$$\quad b=',num2str(b),'\quad c=',num2str(c),'\quad f(x)=x^k, k=',num2str(k),'$$ \quad dFmt=',dFmt,'}'],...
        ['\centerline{cutOff=[',num2str(cutOff(1:N)'),'] \quad $$\sigma=$$',num2str(sigma),'}']},'interpreter','latex','HorizontalAlignment','center');
xlabel('$$x$$','interpreter','latex');
ylabel(ax(1),'$$u(x)$$','interpreter','latex','color','black');
ylabel(ax(2),'Error','color','black');
set(ax(1),'fontsize',12,'xlim',[0,1],'ylim',[0,1.2],'Ycolor','black');
set(ax(2),'fontsize',12,'xlim',[0,1],'Ycolor','black');

