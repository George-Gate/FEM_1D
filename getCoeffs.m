%% Generate the coefficient matrices up to second derivative for boundary-value problem
% Required Input: n, k (f(x)), dFmt, meshType
% Optional Input: extraPointPos, meshWidth
% Output: S - the matrix stands for second derivative of the unknown function
%         C - the matrix stands for first derivative of the unknown function
%         M - the matrix stands for the unknown function itself
%      vecf - the vector stands for the function on the r.h.s. of equation
%     xList - x coordinate of mesh points, excluding boundary points.
%     hList - the width of each mesh grid. 
%                          hList=[xList;1]-[0;xList];
%         N - The total number of mesh grids. 
%                          N=length(xList)+1=length(hList);
%-----------------------------------------------------------------------------------------------------------------------
% Possible dFmt: 
%        FEM - 1D Linear Finite Element Method
%    central - 1D Finite Difference Method, central diference format
%    forward - 1D Finite Difference Method, forward diference format
%   backward - 1D Finite Difference Method, backward diference format
% FEM+Spectrum - 1D Linear Finite Element Method plus Legendre function Spectrum method on each interval
%                Require additional input: cutOff, a N x 1 vector setting the cut off of Spectrum method on each
%                interval.
%                Special Output: fun2id, id2fun, two mapping relation.
%-----------------------------------------------------------------------------------------------------------------------
% Possible meshType: 
%         uniform - (n-1) point uniform mesh on [0,1], i.e. 0=x_0<x_1<...<x_n=1
%       uniformP1 - uniform mesh, plus 1 point on [x_(n-1),x_n]. The position of the extra point can be specified by
%                   variable 'extraPointPos' which is a double between 0 and 1. If extraPointPos doesn't exist, the  
%                   point will be put randomly.
%       uniformP2 - uniform mesh, plus 2 point on two side. The position of the extra point can be specified by
%                   variable 'extraPointPos' which is a 2 x 1 double array between 0 and 1. If extraPointPos doesn't   
%                   exist, the points will be put randomly and symmetrically.
%        shishkin - shishkin mesh where the dense part is near x=1. (n-1) points in each area. The width of dense area
%                   is specified by variable 'meshWidth'.
%   2sideShishkin - two side shishkin mesh that has dense parts near both x=0 and x=1. (n-1) points in each of the
%                   three areas. The width of dense area is specified by variable 'meshWidth'.
%
%

%% generate hList based in meshType and n
clear xList
switch meshType
    case 'uniform'
    % uniform mesh
        hList=ones(n,1)/n;
        xList=linspace(0,1,n+1)';
        xList=xList(2:end-1);
    case 'uniformP1'
    % uniform mesh plus one point on the last mesh grid
        if ~exist('extraPointPos','var') || abs(extraPointPos(1)-0.5)>0.5
            eP=rand();
        else
            eP=extraPointPos(1);
        end
        hList=[ones(n-1,1)/n;(1-eP)/n;eP/n];
        xList=cumsum(hList(1:end-1));
        clear eP;
    case 'uniformP2'
    % uniform mesh plus two points on both side
        if ~exist('extraPointPos','var') || length(extraPointPos)<2 || sum(abs(extraPointPos(1:2)-0.5)>0.5)
            eP1=rand();eP2=eP1;
        else
            eP1=extraPointPos(1);
            eP2=extraPointPos(2);
        end
        hList=[eP1/n;(1-eP1)/n;ones(n-2,1)/n;(1-eP2)/n;eP2/n];
        xList=cumsum(hList(1:end-1));
        clear eP1 eP2;
    case 'shishkin'
    % shishkin mesh
        if ~exist('meshWidth','var') || abs(meshWidth-0.25)>0.25
            error('Invalid meshWidth. meshWidth should be smaller than 1/2.');
        end
        hList=[ones(n,1)*(1-meshWidth)/n;ones(n,1)*meshWidth/n];
        xList=linspace(0,1-meshWidth,n+1)';
        xList=[xList(1:end-1);linspace(1-meshWidth,1,n+1)'];
        xList=xList(2:end-1);
    case '2sideShishkin'
    % two side shishkin mesh
        if ~exist('meshWidth','var') || abs(meshWidth-1/6)>1/6
            error('Invalid meshWidth. meshWidth should be smaller than 1/3.');
        end
        hList=[ones(n,1)*meshWidth/n;ones(n,1)*(1-2*meshWidth)/n;ones(n,1)*meshWidth/n];
        xList=linspace(0,meshWidth,n+1)';
        xList=[xList(1:end-1);linspace(meshWidth,1-meshWidth,n+1)'];
        xList=[xList(1:end-1);linspace(1-meshWidth,1,n+1)'];
        xList=xList(2:end-1);
    otherwise
        error(['Unknow mesh type:',meshType]);
end
if (abs(sum(hList)-1)>1e-10)
    error(['Mesh error: 1-sum(hList)=',num2str(1-sum(hList))]);
end
N=length(hList);
if (N~=length(xList)+1)
    error('Length of xList and hList do not match.');
end
if (xList(end)>1)
    error(['Mesh error: 1-xList(end)=',num2str(1-xList(end))]);
end
% plot([xList;1],hList./hList,'o');

%% generate coefficient matrices based on dFmt and hList
if strcmp(dFmt,'FEM')
% use linear finite element method
%     S=sparse(  diag( ( 1./hList(1:end-1)  +  1./hList(2:end) ).*ones(N-1,1) )  ...
%               +diag(  -1./hList(2:end-1)                      .*ones(N-2,1),-1) ...
%               +diag(  -1./hList(2:end-1)                      .*ones(N-2,1),+1) );
%     C=sparse(  diag(      -1/2                                .*ones(N-2,1),-1) ...
%               +diag(        1/2                                .*ones(N-2,1),+1));
%     M=sparse(  diag( ( hList(1:end-1)  +  hList(2:end) )/3    .*ones(N-1,1) )  ...
%               +diag(   hList(2:end-1)/6                       .*ones(N-2,1),-1) ...
%               +diag(   hList(2:end-1)/6                       .*ones(N-2,1),+1) );
    S= sparse( (1:N-1)',(1:N-1)',( 1./hList(1:end-1)  +  1./hList(2:end) ),N-1,N-1  )  ...
      +sparse( (2:N-1)',(1:N-2)', -1./hList(2:end-1)                      ,N-1,N-1  ) ...
      +sparse( (1:N-2)',(2:N-1)', -1./hList(2:end-1)                      ,N-1,N-1  );
    C= sparse( (2:N-1)',(1:N-2)', -1/2                       ,N-1,N-1  ) ...
      +sparse( (1:N-2)',(2:N-1)',  1/2                      ,N-1,N-1  );
    M= sparse( (1:N-1)',(1:N-1)',( hList(1:end-1)  +  hList(2:end) )/3 ,N-1,N-1  )  ...
      +sparse( (2:N-1)',(1:N-2)',  hList(2:end-1)/6                    ,N-1,N-1  ) ...
      +sparse( (1:N-2)',(2:N-1)',  hList(2:end-1)/6                    ,N-1,N-1  );

%     tic;
    % calc load vector vecf
    vecf=zeros(N-1,1);
    hList2=hList(2:end);
    parfor j=1:N-1
        h1=hList(j);   h2=hList2(j);  x0=xList(j);
        vecf(j)=integral(@(t)(  h1*f( x0+h1*t  )  +  h2*f( x0-h2*t )  ).*(1+t),  -1,  0);
    end
%     disp(['Time used to integrate vecf: ',num2str(toc),'s']);
    % check the symmetry of coeff matrices
    if max(max(abs(S-S')))>eps || max(max(abs(C+C')))>eps || max(max(abs(M-M')))>eps
        error('Symmetry test of the coeff matrices failed.');
    end
    
elseif strcmp(dFmt,'central') || strcmp(dFmt,'forward') || strcmp(dFmt,'backward')
% use finite difference method
    % add a minus sign before to fit with FEM. Precision: O(h1-h2) when h1/=h2
%     S=-sparse(  diag( -( 1./hList(1:end-1)  +  1./hList(2:end) )      .*ones(N-1,1) )  ...
%                +diag(    1./hList(2:end-1)                            .*ones(N-2,1),-1) ...
%                +diag(    1./hList(2:end-1)                            .*ones(N-2,1),+1) );
    S=-sparse( (1:N-1)',(1:N-1)',-( 1./hList(1:end-1)  +  1./hList(2:end) ),N-1,N-1  )  ...
      -sparse( (2:N-1)',(1:N-2)', 1./hList(2:end-1)                        ,N-1,N-1  ) ...
      -sparse( (1:N-2)',(2:N-1)', 1./hList(2:end-1)                        ,N-1,N-1  );
    switch dFmt
        case 'central'
%             C=sparse(  diag( ( hList(2:end)./hList(1:end-1)  -  hList(1:end-1)./hList(2:end)  )/2    .*ones(N-1,1) )  ...
%                       +diag(  -hList(3:end)  ./hList(2:end-1)/2                                      .*ones(N-2,1),-1) ...
%                       +diag(   hList(1:end-2)./hList(2:end-1)/2                                      .*ones(N-2,1),+1) );
            C= sparse( (1:N-1)',(1:N-1)',( hList(2:end)./hList(1:end-1)  -  hList(1:end-1)./hList(2:end)  )/2   ,N-1,N-1  )  ...
              +sparse( (2:N-1)',(1:N-2)',   -hList(3:end)  ./hList(2:end-1)/2                                   ,N-1,N-1  ) ...
              +sparse( (1:N-2)',(2:N-1)',   hList(1:end-2)./hList(2:end-1)/2                                    ,N-1,N-1  );
        case 'forward'
%             C=sparse(  diag( -( hList(1:end-1)./hList(2:end) + 1 )/2     .*ones(N-1,1) )  ...
%                       +diag(  ( hList(1:end-2)./hList(2:end-1) + 1 )/2   .*ones(N-2,1),+1) );
            C= sparse( (1:N-1)',(1:N-1)',  -( hList(1:end-1)./hList(2:end) + 1 )/2    ,N-1,N-1  )  ...
              +sparse( (1:N-2)',(2:N-1)',  ( hList(1:end-2)./hList(2:end-1) + 1 )/2   ,N-1,N-1  );
        case 'backward'
%             C=sparse(  diag(  ( hList(2:end)./hList(1:end-1) + 1 )/2   .*ones(N-1,1) )  ...
%                       +diag( -( hList(3:end)./hList(2:end-1) + 1 )/2   .*ones(N-2,1),-1) );
            C= sparse( (1:N-1)',(1:N-1)',    ( hList(2:end)./hList(1:end-1) + 1 )/2     ,N-1,N-1  )  ...
              +sparse( (2:N-1)',(1:N-2)',   -( hList(3:end)./hList(2:end-1) + 1 )/2     ,N-1,N-1  );
    end
%     M=sparse(  diag(  (hList(1:end-1) + hList(2:end))/2   .*ones(N-1,1) )  );
    M=sparse( (1:N-1)',(1:N-1)', (hList(1:end-1) + hList(2:end))/2, N-1, N-1);
    

    %calc f
    vecf=f(xList).*(hList(1:end-1)+hList(2:end))/2;
elseif strcmp(dFmt,'FEM+Spectrum')
    if ~exist('cutOff','var') || length(cutOff)<N
        error('Invalid argument: cutOff.');
    end
    tic;
    % mapping of base function and vector index
    fun2id.phi=(1:N-1)';   % linear base function of FEM mesh
    fun2id.psi=cell(N,1);  % fun2id.psi{k} is the Lobatto Function on interval [x_(k-1),x_k], fun2id.psi{k}(j)-->Lo_(j+1)
    counter=N-1;
    for m=1:N
        fun2id.psi{m}=counter+(1:cutOff(m))';
        counter=counter+cutOff(m);
    end
    
    id2fun=zeros(counter,2);  % id2fun(1:N-1,1) the index of linear base function, id2fun(1:N-1,2) --> 0
                              % id2fun(N:end,1) the interval index of Lobatto Function
                              % id2fun(N:end,2) the order of Lobatto Function
    id2fun(1:N-1,1)=(1:N-1);
    counter=N-1;
    for m=1:N
        id2fun(counter+(1:cutOff(m))',1)=m;
        id2fun(counter+(1:cutOff(m))',2)=(1:cutOff(m))';
        counter=counter+cutOff(m);
    end
    
    clear counter;
    
    % diag blocks
    % S00, C00 and M00
    S=sparse(  diag( ( 1./hList(1:end-1)  +  1./hList(2:end) ), 0) ...  .*ones(N-1,1)
              +diag(  -1./hList(2:end-1)                      ,-1) ...  .*ones(N-2,1)
              +diag(  -1./hList(2:end-1)                      ,+1) ); % .*ones(N-2,1)
    C=sparse(  diag(      -1/2                                .*ones(N-2,1),-1) ...
              +diag(       1/2                                .*ones(N-2,1),+1));
    M=sparse(  diag( ( hList(1:end-1)  +  hList(2:end) )/3    , 0) ...  .*ones(N-1,1)
              +diag(   hList(2:end-1)/6                       ,-1) ...  .*ones(N-2,1)
              +diag(   hList(2:end-1)/6                       ,+1) ); % .*ones(N-2,1)
          
    max(abs(S(:)-S1(:)))
    max(abs(C(:)-C1(:)))
    max(abs(M(:)-M1(:)))
    
          
    % Smm, Cmm and Mmm
    ub=max(cutOff(1:N));   % max upper bound of the order of legendre polynomial
    id=(1:ub)';
%     Smmhm=sparse(  diag( 2*ones(ub  ,1) ));
    Smmhm=sparse((1:ub)',(1:ub)',2,ub,ub);
%     Cmm2  =sparse(  diag( -1./sqrt((2*id(2:end  )+1).*(2*id(2:end  )-1)) ,-1) ... .*ones(ub-1,1)
%                 +diag(  1./sqrt((2*id(1:end-1)+1).*(2*id(1:end-1)+3)) ,+1));    %.*ones(ub-1,1)
    Cmm  =sparse((2:ub)',(1:ub-1)',-1./sqrt((2*id(2:end  )+1).*(2*id(2:end  )-1)),ub,ub) ... 
         +sparse((1:ub-1)',(2:ub)',  1./sqrt((2*id(1:end-1)+1).*(2*id(1:end-1)+3)),ub,ub); 
%     Mmm_hm=0.5*sparse(  diag((  (4*id+2)./( (2*id-1).*(2*id+1).*(2*id+3) )                             ) )  ...  .*ones(ub  ,1)
%                        +diag((  -1./( (2*id(3:end  )-1).*sqrt((2*id(3:end  )-3).*(2*id(3:end  )+1)) )  ),-2) ... .*ones(ub-2,1)
%                        +diag((  -1./( (2*id(1:end-2)+3).*sqrt((2*id(1:end-2)+1).*(2*id(1:end-2)+5)) )  ),+2) ); %.*ones(ub-2,1)
    Mmm_hm=0.5*(sparse((1:ub)',(1:ub)',   (4*id+2)./( (2*id-1).*(2*id+1).*(2*id+3) )                           ,ub,ub  )  ... 
               +sparse((3:ub)',(1:ub-2)', -1./( (2*id(3:end  )-1).*sqrt((2*id(3:end  )-3).*(2*id(3:end  )+1)) ),ub,ub  ) ... 
               +sparse((1:ub-2)',(3:ub)', -1./( (2*id(1:end-2)+3).*sqrt((2*id(1:end-2)+1).*(2*id(1:end-2)+5)) ),ub,ub  ) ); 

    for m=1:N
        ub=cutOff(m);   % upper bound of the order of legendre polynomial
        S=blkdiag(S,Smmhm(1:ub,1:ub)/hList(m));
        C=blkdiag(C,Cmm(1:ub,1:ub));
        M=blkdiag(M,Mmm_hm(1:ub,1:ub)*hList(m));
    end
    
    % off diag blocks
    tmp1=1/sqrt(6);
    tmp2=sqrt(2/3)/4;
    tmp3=sqrt(2/5)/12;
    for m=1:N
        k1=fun2id.psi{m}(1);
        k2=fun2id.psi{m}(2);
        if (m<N)
            C(m,k1)=tmp1;                    C(k1,m)=-C(m,k1);
            M(m,k1)=-hList(m)*tmp2;          M(k1,m)=M(m,k1);
            M(m,k2)=-hList(m)*tmp3;          M(k2,m)=M(m,k2);
        end
        if (m>1)
            C(m-1,k1)=-tmp1;                 C(k1,m-1)=-C(m-1,k1);
            M(m-1,k1)=-hList(m)*tmp2;        M(k1,m-1)=M(m-1,k1);
            M(m-1,k2)=+hList(m)*tmp3;        M(k2,m-1)=M(m-1,k2);
        end
    end
    
%     disp(['Time used to assemble coefficient matrices: ',num2str(toc),'s']);
    % check the symmetry of coeff matrices
    if max(max(abs(S-S')))>eps || max(max(abs(C+C')))>eps || max(max(abs(M-M')))>eps
        error('Symmetry test of the coeff matrices failed.');
    end
    tic;
    % calc load vector
    vecf=zeros(N-1+sum(cutOff(1:N)),1);
    % for phi
    for j=1:N-1
        h1=hList(j);   h2=hList(j+1);  x0=xList(j);
        vecf(j)=integral(@(t)(  h1*f( x0+h1*t  )  +  h2*f( x0-h2*t )  ).*(1+t),  -1,  0);
    end
%     disp(['Time used to integrate f(1:N-1): ',num2str(toc),'s']);
    tic;
    % for psi (Lobatto Function)
    % use the order information of f(x) ,i.e. k, to accelerate integration.
    for m=1:N
        for order=1:min(k+1,cutOff(m))
            hm_2=hList(m)/2;
            if (m>1)
                xm_1=xList(m-1);
            else
                xm_1=0;
            end
            vecf(fun2id.psi{m}(order))=hm_2/( sqrt(4*order+2) )*integral(@(t)f(hm_2*(t+1)+xm_1).*(legendreP_N(order+1,t)-legendreP_N(order-1,t)),-1,1);  
        end
    end
%     disp(['Time used to integrate f(N:end): ',num2str(toc),'s']);
    clear Smmhm Cmm Mmm_hm id ub m tmp1 tmp2 tmp3 hm_2 xm_1 order
else
    error(['Unknow differential format:',dFmt]);
end