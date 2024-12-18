function [zero,res,niter]=multinewton(f,df,x0,tol,nmax,silent,varargin)
% MULTINEWTON Find function zeros using Newtons method in more dimensions.
%
%   ZERO=MULTINEWTON(FUN,DFUN,X0,TOL,NMAX,SILENT) 
%   tries to find the zero ZERO of the continuous and 
%   differentiable function FUN nearest to X0 using Newtons method. 
%   FUN is a function from R^n to R^n (a vector valued function)
%   Its derivative is DFUN, which a function from
%   R^n to R^(n^2), i.e., a matrix valued function.
%   Both accept a column vector input (and FUN also returns a column vector).
%   FUN and DFUN can also be inline objects or anonymous functions.
%
%   If the search fails an error message is displayed.
%   If SILENT=1 then no error messages are displayed (SILENT=0 by default).
%
%   The tolerance TOL and the maximum number NMAX of iterations
%   have default values 1e-14 and 100.
%   Convergence test is: norm(x(n)-x(n-1)) < TOL.
%
%   [ZERO,RES,NITER]= MULTINEWTON(FUN,...) 
%   returns the value of the residual RES (as a vector)
%   in ZERO and the number of itererations NITER.

% initialization
if ~exist('nmax')
  nmax=1000;
end
if ~exist('tol')
  tol=1e-14;
end
if ~exist('silent')
  silent=0;  
end

x = x0; 
niter = 0; 
diff = tol+1;
while norm(diff) >= tol & niter <= nmax
   fx = feval(f,x,varargin{:}); 
   dfx = feval(df,x,varargin{:});
   niter = niter + 1;
   diff = - dfx\fx;
   norm(diff);
   x = x + diff;     
end
if niter > nmax & silent==0
    fprintf(['newton stopped without converging to the desired tolerance\n',...
    'because the maximum number of iterations was reached\n']);
end
zero = x;
res = feval(f,x,varargin{:});

return
