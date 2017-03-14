%% Generate the coefficient matrices
% Required Input: n, k (f(x)), dFmt, meshType
% Optional Input: extraPointPos, meshWidth
% Output: S - the matrix stands for second derivative of the unknown function
%         C - the matrix stands for first derivative of the unknown function
%         M - the matrix stands for the unknown function itself
%      vecf - the vector stands for the function on the r.h.s. of equation
%     xList - x coordinate of mesh points, excluding boundary points.
%     hList - the width of each mesh grid. 
%                          hList=[xList;1]-[0;xList];
%         N - The total number of mesh points, excluding boundary points. 
%                          N=length(xList);
% Possible dFmt: 
%        FEM - Finite Element Method
%    central - Finite Difference Method, central diference format
%    forward - Finite Difference Method, forward diference format
%   backward - Finite Difference Method, backward diference format
% Possible meshType: 
%         uniform - (n-1) point uniform mesh on [0,1], i.e. 0=x_0<x_1<...<x_n=1
%       uniformP1 - uniform mesh, plus 1 point on [x_(n-1),x_n]. The position of the extra point can be specified by
%                   variable 'extraPointPos' which is a double between 0 and 1. If extraPointPos doesn't exist, the  
%                   point will be put randomly.
%        shishkin - shishkin mesh where the dense part is near x=1. (n-1) points in each area. The width of dense area
%                   is specified by variable 'meshWidth'.
%   2sideShishkin - two side shishkin mesh that has dense parts near both x=0 and x=1. (n-1) points in each of the
%                   three areas. The width of dense area is specified by variable 'meshWidth'.
%
%

%% generate hList based in meshType and n
switch meshType
    case 'uniform'
    % uniform mesh
        hList=ones(n,1)/n;
    case 'uniformP1'
    % uniform mesh plus one point on the last mesh grid
        if ~exist('extraPointPos','var') || abs(extraPointPos-0.5)>0.5
            extraPointPos=rand();
        end
        hList=[ones(n-1,1)/n;extraPointPos/n;(1-extraPointPos)/n];
    case 'shishkin'
    % shishkin mesh
        if ~exist('meshWidth','var') || abs(meshWidth-0.25)>0.25
            error('Invalid meshWidth. meshWidth should be smaller than 1/2.');
        end
        hList=[ones(n,1)*(1-meshWidth)/n;ones(n,1)*meshWidth/n];
    case '2sideShishkin'
    % two side shishkin mesh
        if ~exist('meshWidth','var') || abs(meshWidth-1/3)>1/6
            error('Invalid meshWidth. meshWidth should be smaller than 1/3.');
        end
        hList=[ones(n,1)*meshWidth/n;ones(n,1)*(1-2*meshWidth)/n;ones(n,1)*meshWidth/n];
    otherwise
        error(['Unknow mesh type:',meshType]);
end
if (abs(sum(hList)-1)>1e-10)
    error(['Mesh error: 1-sum(hList)=',num2str(1-sum(hList))]);
end
xList=cumsum(hList(1:end-1));
N=length(xList);
% plot([xList;1],hList./hList,'o');

%% generate coefficient matrices based on dFmt and hList
if strcmp(dFmt,'FEM')
% use linear finite element method
    S=sparse(  diag( ( 1./hList(1:end-1)  +  1./hList(2:end) ).*ones(N  ,1) )  ...
              +diag(  -1./hList(2:end-1)                        .*ones(N-1,1),-1) ...
              +diag(  -1./hList(2:end-1)                      .*ones(N-1,1),+1) );
    C=sparse(  diag(      -1/2                                .*ones(N-1,1),-1) ...
              +diag(       1/2                                .*ones(N-1,1),+1));
    M=sparse(  diag( ( hList(1:end-1)  +  hList(2:end) )/3    .*ones(N  ,1) )  ...
              +diag(   hList(2:end-1)/6                         .*ones(N-1,1),-1) ...
              +diag(   hList(2:end-1)/6                       .*ones(N-1,1),+1) );
    % calc f
    vecf=zeros(N,1);
    for j=1:N
        h1=hList(j);   h2=hList(j+1);  x0=xList(j);
        vecf(j)=integral(@(t)(  h1*f( x0+h1*t  )  +  h2*f( x0-h2*t )  ).*(1+t),  -1,  0);
    end
    
elseif strcmp(dFmt,'central') || strcmp(dFmt,'forward') || strcmp(dFmt,'backward')
% use finite difference method
    % add a minus sign before to fit with FEM. Precision: O(h1-h2) when h1/=h2
    S=-sparse(  diag( -( 1./hList(1:end-1)  +  1./hList(2:end) )      .*ones(N  ,1) )  ...
               +diag(    1./hList(2:end-1)                            .*ones(N-1,1),-1) ...
               +diag(    1./hList(2:end-1)                            .*ones(N-1,1),+1) );
    switch dFmt
        case 'central'
            C=sparse(  diag( ( hList(2:end)./hList(1:end-1)  -  hList(1:end-1)./hList(2:end) )/2   .*ones(N  ,1) )  ...
                      +diag(  -hList(3:end)  ./hList(2:end-1)/2                                    .*ones(N-1,1),-1) ...
                      +diag(   hList(1:end-2)./hList(2:end-1)/2                                    .*ones(N-1,1),+1) );
        case 'forward'
            C=sparse(  diag( -( hList(1:end-1)./hList(2:end) + 1 )/2     .*ones(N  ,1) )  ...
                      +diag(  ( hList(1:end-2)./hList(2:end-1) + 1 )/2   .*ones(N-1,1),+1) );
        case 'backward'
            C=sparse(  diag(  ( hList(2:end)./hList(1:end-1) + 1 )/2   .*ones(N  ,1) )  ...
                      +diag( -( hList(3:end)./hList(2:end-1) + 1 )/2   .*ones(N-1,1),-1) );
    end
    M=sparse(  diag(  (hList(1:end-1) + hList(2:end))/2   .*ones(N,1) )  );

    %calc f
    vecf=f(xList).*(hList(1:end-1)+hList(2:end))/2;
else
    error(['Unknow differential format:',dFmt]);
end