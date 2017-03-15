function [ x , y ] = linkIntersection( xList, line1, line2 )
%Find the intersection of two line.
%% Detailed explanation
%  Input: xList, line1 and line2, three column vector with same dimenesion. 
%  Output: x and y, the x, y coordinate of intersection line.
%
% <<..\documentation\linkIntersection\how_to_find_point.jpg>>
% 
%% Test code
%   xList=(0:0.1:1)';
%   line1=[0;2;2;0;3;1;2;3;2;1;3];
%   line2=[0;1;1;2;1;2;3;2;1;0;3];
%   [intX,intY]=linkIntersection(xList,line1,line2);
%   plot(xList,line1,'-o',xList,line2,'-o',intX,intY,'--*');


%% Code
    if length(xList)~=length(line1) || length(line1)~=length(line2)
        error('xList, line1 and line2 should be three vectors with same length.');
    end
    len=length(xList);
    x=zeros(len-1,1);
    y=x;
    topX=1;
    
    sgnList=sign(line1-line2)+0.5;  % treat line1==line2 as line1>line2
    for i=1:len-1
        if ( sgnList(i+1)*sgnList(i)<0 )
            % there is an intersection between point i and i+1
            d=line1(i:i+1)-line2(i:i+1);
            x0=( d(1)*xList(i+1) - d(2)*xList(i) )/( d(1)-d(2) );
            x(topX)=x0;
            y(topX)=(  line1(i+1)*( x0-xList(i) )  +  line1(i)*( xList(i+1)-x0 )   )/( xList(i+1)-xList(i) );
            topX=topX+1;
        else
            % there is no intersection between point i and i+1
            % check if there is an intersection point on the neighbour interval
            if i==1
                % for the first point
                x(topX)=xList(1);
                y(topX)=(line1(1)+line2(1))/2;
                topX=topX+1;
            else
                % for inner points, if the fromer interval also have no intersection, then ...
                if (sgnList(i)*sgnList(i-1)>0)
                    x(topX)=xList(i);
                    y(topX)=(line1(i)+line2(i))/2;
                    topX=topX+1;
                end
            end
            if i==len-1
                % for the last point
                x(topX)=xList(len);
                y(topX)=(line1(len)+line2(len))/2;
                topX=topX+1;
            end
        end
    end
    % trim x and y
    x=x(1:topX-1);
    y=y(1:topX-1);

end



