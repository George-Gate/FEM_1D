%% analytical solution
% Required Input: epsilon, b, c, k (f(x))
% Output: @(x)anaSol(x)

if (abs(c)>eps)
    u0={@(x)1/c,...
        @(x)x/c-b/c^2,...
        @(x)x.^2/c-2*b*x/c^2+2*epsilon/c^2+2*b^2/c^3,...
        @(x)x.^3/c-3*b*x.^2/c^2+6/c^2*(epsilon+b^2/c)*x-6*b/c^3*(2*epsilon+b^2/c)};
    lambda_p=( b+sqrt(b^2+4*epsilon*c) )/2/epsilon;
    lambda_m=( b-sqrt(b^2+4*epsilon*c) )/2/epsilon;
else
    e2b=epsilon/b;
    u0={@(x)x/b,...
        @(x)  (x/2/b+e2b/b).*x,...
        @(x) ((x/3/b+e2b/b).*x+2*(e2b)^2/b).*x,...
        @(x)(((x/4/b+e2b/b).*x+3*(e2b)^2/b).*x+6*(e2b)^3/b).*x};
    lambda_p=( b+sqrt(b^2) )/2/epsilon;
    lambda_m=( b-sqrt(b^2) )/2/epsilon;
end
    

% solve for coefficients
if (lambda_p-lambda_m<300)
    xx=exp(lambda_p);yy=exp(lambda_m);
    logCoeff=log( [ yy*u0{k+1}(0)-u0{k+1}(1) ; u0{k+1}(1)-xx*u0{k+1}(0) ] ) - log(xx-yy);
else
    % for the case where exp(lambda_p)/exp(lambda_m)>>1
    yy=exp(lambda_m);
    logCoeff=[log( yy*u0{k+1}(0) - u0{k+1}(1))-lambda_p ; log(-u0{k+1}(0))];
end

anaSol=@(x)exp(logCoeff(1)+lambda_p*x)+exp(logCoeff(2)+lambda_m*x)+u0{k+1}(x);
clear xx yy u0 logCoeff lambda_p lambda_m e2b