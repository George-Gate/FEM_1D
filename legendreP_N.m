function [ result ] = legendreP_N( nList,xList )
%% Calculate the value of Legendre Polynomials at given position
%  [Usage]
%       result=legendreP_N( n, x )
%   n and x should be vectors and n>=0
%       if x is a col vector,
%             result=[L( n(1) , x ),L( n(2) , x ),L( n(3) , x ), ... ,L( n(end) , x )]
%       if x is a row vector,
%             result=[L( n(1) , x );L( n(2) , x );L( n(3) , x ); ... ;L( n(end) , x )]
%

    if ~isvector(xList) || ~isvector(nList)
        error('nList and xList must be vectors.');
    end
    maxN=max(nList);
    P=cell(maxN+1,1);
    P{1}=ones(size(xList));
    P{2}=xList;
    for n=1:maxN-1
        P{n+2}=( (2*n+1)*xList.*P{n+1}-n*P{n} )/(n+1);
    end
    % set output
    if (size(xList,1)==1)
        result=zeros(length(nList),length(xList));
        for i=1:length(nList)
            result(i,:)=P{nList(i)+1};
        end
    elseif (size(xList,2)==1)
        result=zeros(length(xList),length(nList));
        for i=1:length(nList)
            result(:,i)=P{nList(i)+1};
        end
    end
end

