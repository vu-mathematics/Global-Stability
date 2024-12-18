function [t,x,v,J,diffxopt] = daes_CN(odedae,minnewton,compute_functions,Tmax);

% function [t,x,x0,Y0,Finv] = daes_CN(odedae,minnewton,compute_functions,Tmax);
%
% Belongs to "Maintaining maximal metabolic flux using gene expression control" 
% by Planque et al. 
% 
% This function integrates the metabolic network shown below, and corresponds 
% to Figure 2 in the paper. It implements the DyRAC framework for one specific pathway. 
%
% Arguments are:
% odedae: solve as ODE or as DAE.  odedae==1 (ODE) is default. Set to 0 otherwise
% minnewton: find start optimum using newton's method (default) or minimization (set to zero)
% compute_functions: (re-)compute the symbolic functions, Jacobian and 'Hessian', which 
% 	are stored in a separate directory after computation, for easy reuse next time.
% Tmax: the ODE/DAE is integrated for three subsequent intervals of length Tmax.
%
% x1 -- x2 -\
%            \
%			  -- x4 ---> OUT
%		     /
% N --- x3 -/
% 
% means
% 
% C ---- C'-\
%            \
%			  -- CN ---> output
%		     /
% N ---- N'-/
%
% N is treated as parameter, not as dynamic variable in this example

if nargin < 4
  Tmax = 50;
elseif nargin < 3
  compute_functions=0; 
elseif nargin < 2
  minnewton = 1; 
elseif nargin < 1
  odedae=1;  
end

t = [];x = [];
format compact
n = 4; % number of metabolites. x5 is not a dynamic variable. it is a parameter here
nenz = 4; % number of reactions/enzymes
% parameters used in the kinetics functions
A = [
3 1 3 3 1 2 3 1 3 2; 
5 1 3 3 2 3 2 1 3 4 ; 
4 6 5 1 3 4 1 3 3 1;
3 2 1 4 1 2 3 2 1 2;
4 1 6 2 1 5 6 2 3 9;
1 3 1 2 1 43 12 1 2 2;
1 4 1 1 4 2 1 2 3 4];
BCs = [5]; %externally fixed concentrations x1 and x6. x2 is fixed but as parameter
sensors = [2]; % defines which metabolite is the sensor 

% stitching parts together
fullxvector=@(x,BCs) [BCs(1) x(1:end)'];

% the stoichiometry matrix
N = [0  0  0  0; %x1 is not dynamic
	 1  0 -1  0;
	 0  1 -1  0;
	 0  0  1 -1];

e = [0.23 0.3 0.1 0.2];

es = sum(e); e=e/es; % normalise to sum=1
NN = 4; % external N concentration
OUT = 0.01;
params{1} = A;
params{2} = e;
params{3} = n;
params{4} = BCs;
params{5} = sensors;
params{6} = N;
params{7} = nenz;
params{8} = odedae;
params{9} = NN;

x0 = 10*rand(n-length(BCs),1);
x0 = [5 3 4]';
X0 = fullxvector(x0,BCs);


if(compute_functions) 

  %define the object function symbolically, and compute its symbolic Jacobian, 
  % which yields the optimum equations.
  syms x1 x2 x3 x4;
  vars=[x1 x2 x3 x4];
  f1=symfun((A(1,1)*x1/A(2,1)-A(3,1)*x2/A(4,1))/...			% facilitated diffusion kinetics
	((1 + x1/A(2,1))*(1 + x2/A(4,1))), vars);
  f2=symfun((A(1,2)*NN/A(2,2)-A(3,2)*x3/A(4,2))/...			% facilitated diffusion kinetics
	((1 + NN/A(2,2))*(1 + x3/A(4,2))), vars);
  f3 = symfun((A(1,3)*x2*x3/A(2,3)/A(3,3) - A(4,3)*x4/A(5,3))/...
	((1 + x2/A(2,3) + x4/A(5,3))*(1 + x3/A(3,3))) ,vars);  
  f4=symfun((A(1,4)*x4/A(2,4)- A(3,4)*OUT/A(4,4))/(1 + x4/A(2,4) + OUT/A(4,4)), vars);

  % compute the Jacobians of the various objects
  disp('... computing the Jacobians....')
  FTNS = symfun([f1;f2;f3;f4],vars);	% the various kinetic functions
  OBJ_PARTS = symfun([1/f1; 1/f2; 1/f3; 1/f4],vars);	% the parts of the objective function
  OBJSYM = symfun(1/f1 + 1/f2 + 1/f3 + 1/f4, vars);		% the objective function
  DOBJDX = symfun(jacobian(OBJSYM,[x2 x3 x4]),vars); 	% the equations in optimum
  DOBJDXX = symfun(jacobian(DOBJDX,vars),vars);			% used in the explicit ODE description, when the DAE has been converted to an ODE.
  DOBJDXX_newton = symfun(jacobian(DOBJDX,[x2 x3 x4]),vars); % used to find an initial solution to the DAE system, using Newton's method

  % the symbolic functions are converted to matlab functions for fast execution
  ftns1      = matlabFunction(FTNS);
  ftns       = @(x) ftns1(x(1),x(2),x(3),x(4));
  Obj_parts1 = matlabFunction(OBJ_PARTS);
  Obj_parts  = @(x) Obj_parts1(x(1),x(2),x(3),x(4));
  Objsym1    = matlabFunction(OBJSYM);
  Objsym     = @(x) Objsym1(x(1),x(2),x(3),x(4));
  dObjdx1    = matlabFunction(DOBJDX);
  dObjdx     = @(x) dObjdx1(x(1),x(2),x(3),x(4));
  dObjdxx1   = matlabFunction(DOBJDXX);
  dObjdxx     = @(x) dObjdxx1(x(1),x(2),x(3),x(4));
  dObjdxx_newton1 = matlabFunction(DOBJDXX_newton);
  dObjdxx_newton  = @(x) dObjdxx_newton1(x(1),x(2),x(3),x(4));

  save('CN_1.mat','vars','ftns','Obj_parts','Objsym','dObjdx',...
	'dObjdxx','dObjdxx_newton','A','NN','n','nenz','OUT','sensors');
else
  load('CN_1.mat');
end


disp('... computing the QSSA....')
% first we do a QSSA on x. No enzyme dynamics. The result is a well-ordered 
% concentration vector
Metsys = @(t,x) metsys(t,x,ftns,params);
options = odeset('RelTol',1e-6,'AbsTol',1e-8);
[t,x0] = ode45(Metsys,[0,1000],X0,options);
x0 = x0(end,[2:end])'

disp('... computing the initial optimum....')
%compute the initial optimum

BCs0=[5];
z0=x0;
zopt = compute_optimum(BCs0,z0,fullxvector,Objsym,dObjdx,dObjdxx_newton,...
				A,NN,OUT,minnewton);
x0(sensors-1) = zopt(sensors);

% compute the real optimum
zopt1 = compute_optimum(4,z0,fullxvector,Objsym,dObjdx,dObjdxx_newton,...
				A,NN,OUT,minnewton);
Eopt1 = Obj_parts(zopt1); som=sum(Eopt1);
Fopt1 = ftns(zopt1);
Jopt1 = Eopt1(end)*Fopt1(end) / som;

X0 = fullxvector(x0,BCs);
x0 = [X0(:); e(:); zopt(:)]; % length 3n-1;
x0 = x0(:)'
%x0(3:4)=0;
disp('... integrating the ODE ....')
options = odeset;
if odedae == 0
  % treat the system as a DAE, so define a Mass matrix
  % the difference in the righthand side of the DAE/ODE is defined in the ode() function below
  Mass = zeros(2*n+nenz); % n metabolites's, n-2 enzymes's, n predicted optima elements
  Mass(1:n+nenz,1:n+nenz) = eye(n+nenz);
  options = odeset(options,'Mass',Mass,'MStateDependence','none');
end
options =  odeset(options,'RelTol',1e-10,'AbsTol',1e-10);
ODE = @(t,x) ode(t,x,ftns,Obj_parts,dObjdxx,dObjdx,params);
[t,x] = ode15s(ODE,[0,Tmax],x0,options);
ODE(t(end),x(end,:)')'
J = Jopt1*ones(length(t),1);

% now change the outside concentration, and integrate again
BCs=[2];
params{4}=BCs;
x0=x(end,:);
x0(1)=BCs(1); 
[t1,x1] = ode15s(ODE,[t(end),t(end)+Tmax],x0,options);
t=[t;t1]; x=[x;x1];

% compute the real optimum
zopt2 = x1(end,1:n);
Eopt2 = Obj_parts(zopt2); som=sum(Eopt2);
Fopt2 = ftns(zopt2);
Jopt2 = Eopt2(end)*Fopt2(end) / som
J = [J; Jopt2*ones(length(t1),1)];

% now change the outside concentration once more, and integrate again
BCs=[6];
params{4}=BCs;
x0=x(end,:);
x0(1)=BCs(1); 
[t1,x1] = ode15s(ODE,[t(end),t(end)+Tmax],x0,options);
t=[t;t1]; x=[x;x1];

% compute the real optimum
zopt3 = x1(end,1:n);
Eopt3 = Obj_parts(zopt3); som=sum(Eopt3);
Fopt3 = ftns(zopt3);
Jopt3 = Eopt3(end)*Fopt3(end) / som
J = [J; Jopt3*ones(length(t1),1)];

% compute the flux values in all the time points
for i=1:length(t)
  f = ftns(x(i,1:n));
  v(i,1:nenz) =x(i,n+1:n+nenz).*f';
end

% a plot command which reproduces the figure in the paper
diffxopt = plot_optima_branched(t,x,v,0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dx = metsys(t,x,ftns,params);

% defines the metabolic pathway, for fixed enzyme concentrations

A = params{1};
e = params{2};
n = params{3};
BCs = params{4};
sensors = params{5};
N = params{6};
nenz = params{7};
NN = params{9};

V = e(:).*ftns(x);

dx = N*V(:); 
[t dx' x(:)'];
dx = dx(:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dxe = ode(t,X,ftns,Obj_parts,dObjdxx,dObjdx,params);

% defines the RHS of the ODE/DAE system

A = params{1};
n = params{3};
BCs = params{4};
sensors = params{5};
N = params{6};
nenz = params{7};
odedae = params{8};
nsens = length(sensors);
CC = params{9};

x = X(1:n); % metabolites
e = X(n+1:n+nenz); % enzymzes
z = X(n+nenz+1:end); % raw sensor optimum

% the metabolic pathway: x' = Nv(x) = N e f(x)
V = e(:).*ftns(x);
dx = N*V(:);

% now define RHS of dz/dt (in the DAE case, the LHS is set to zero of course)
dz = choose_dae_ode(odedae,nsens,sensors,dx,n,x,z,dObjdxx,dObjdx);

% define the de/dt RHS's
Eopt = Obj_parts(z); som=sum(Eopt);
for i=1:nenz
  de(i) = 1/10*(Eopt(i)/som - e(i));
end

% glue everything together, and make into a column vector
dxe = [dx(:);de(:);dz(:)];

%%%%%%%%%%%%

function diffxopt = plot_optima_branched(t,y,v,doprint)

% produces the plots in the paper

if nargin == 2
  doprint = 0;
end

n = 4; 
nenz = 4;
x = y(:,1:n);
e = y(:,n+1:n+nenz);
opt = y(:,n+nenz+1:end);
BCs = 1;
rest = setdiff([1:n],BCs);
sensor = 2;
figure(1)
  clf
  hold on
  for i=1:nenz
    p(i) = plot(t,e(:,i),'LineWidth',2);
  end
  ylim = get(gca,'YLim');
  plot([50,50],[0,ylim(2)],'k--')
  plot([100,100],[0,ylim(2)],'k--')
  hold off  
  legend([p],'Enzyme 1', 'Enzyme 2','Enzyme 3','Enzyme 4')
  set(gca,'Box','on') 
%  xlabel('time');
%  ylabel('enzyme concentration');


%  grid on
%  title ('enzymes')
figure(2)
  clf
  hold on
  p2=plot(t,opt,'b-','LineWidth',2);
  p1=plot(t,x(:,rest(2:end)),'r-','LineWidth',2);
  plot(t,x(:,rest(1)),'r--','LineWidth',2);
  p3=plot(t,x(:,BCs),'Color',[0 0.6 0],'LineWidth',2);
  grid off
  set(gca,'YLim',[0,max(y(:))*1.1]);
  ylim = get(gca,'YLim');
  plot([50,50],[0,ylim(2)],'k--')
  plot([100,100],[0,ylim(2)],'k--')
  set(gca,'Box','on') 
  hold off  
%  xlabel('time');
%  ylabel('concentration');
  legend([p1(1) p3(1) p2(1)],'Internal metabolites', 'External Metabolites', 'Predicted optimum')

figure(3)
  clf
  diffxopt = sum((x-opt).^2,2);
  optend = repmat(opt(end,:),length(t),1);
  diffxoptend=sum((x-optend).^2,2);
  hold on
  plot(t,diffxopt,'b-','LineWidth',2);%,t,diffxoptend)
  ylim = get(gca,'YLim');
  plot([50,50],[0,ylim(2)],'k--')
  plot([100,100],[0,ylim(2)],'k--')
  grid off  
  set(gca,'Box','on') 
  hold off
%  xlabel('time');
%  ylabel('| metabolites - predicted optimum |^2');

figure(4)
  clf
  hold on
  plot([0,t(end)],[0,0],'k-')
  for i=1:size(v,2)
    q(i) = plot(t,v(:,i),'LineWidth',2);
  end
  set(gca,'YLim',[min(v(:)*1.1),max(v(:))*1.1]);
  ylim = get(gca,'YLim');
  plot([50,50],[0,ylim(2)],'k--')
  plot([100,100],[0,ylim(2)],'k--')
  set(gca,'Box','on') 
  hold off
  legend([q],'v_1','v_2','v_3','v_4')
%  xlabel('time')
%  ylabel('reaction flux')

if doprint
  print_fig(1, 'pics/CN_1_enz.jpg', [0 0 5 3]);
  print_fig(2, 'pics/CN_1_metab.jpg', [0 0 5 3]);
  print_fig(3, 'pics/CN_1_conv.jpg', [0 0 5 3]);
  print_fig(4, 'pics/CN_1_fluxes.jpg', [0 0 5 3]);

end


%%%%%%%%%%%%%%%%%%

function dz = choose_dae_ode(odedae,nsens,sensors,dx,n,x,z,dObjdxx,dObjdx);

% produces the RHS of the predicted sensor equations
if odedae == 1 % ODE
  % we make the DAE explicit, by converting it to an ODE
  dxsensors = [zeros(n-nsens,1); dx(sensors)];
  DZ = dObjdxx(z);
  [DZr,DZc] = size(DZ);
  F = [DZ;zeros(nsens,n)];
  for i=1:nsens,
    F(DZr+i,sensors(i))=-1; 
  end
  condest(F);
  dz =  -F \ dxsensors;
else % DAE
  % dO/dxi = 0, and xi_S = x_S
  dxy = x(sensors) - z(sensors);
  dz = dObjdx(z);
  dz = [dxy(:); dz(:)];
end



%%%%%%%%%%%%%%%%%%%%%

function zopt=compute_optimum(BCs0,z0,fullxvector,Objsym,dObjdx,dObjdxx_newton,...
				A,NN,OUT,minnewton);

% compute the initial optimum, to start the DAE/ODE system. 

f = @(x) dObjdx(fullxvector(x,BCs0))';
df = @(x) dObjdxx_newton(fullxvector(x,BCs0));
if minnewton == 1
  % just use Newton's method to find a critical point of the objective function
  zopt = multinewton(f,df,z0,1e-7);
else
  % actually minimise the objective function. 
  % this is most fruitfully done in "logarithmic variables", as shown in the paper,
  % since then the objective function is strictly convex. 
  Nconstr = [ 1  0  0;
			 -1  1  0;
			 -1 -1  1;
			  0  0 -1];
  lA = log(A);
  logKeq = [[1 -1 -1 1]*lA(1:4,1) [1 -1 -1 1]*lA(1:4,2) [1 -1 -1 -1 1]*lA(1:5,3) ...
			[1 -1 -1 1]*lA(1:4,4)] 
  Kconstr = logKeq + [log(BCs0(1)),log(NN),0,-log(OUT)]
  % in log-concentration vars the inequality constraints read N' * xi <= log(K);

  % change variables
  xi = log(z0);

  % perform the optimisation
  obj = @(xi) Objsym(fullxvector(exp(xi),BCs0));
  opts = optimoptions('fmincon','TolCon',1e-10,'TolFun',1e-10);
  zopt = fmincon(obj,log(z0),Nconstr,Kconstr,[],[],[],[],[],opts);
  
  % you can check here if the optimum is actually inside the domain by taking
  % away the %'s in the next two lines
  %disp('these constraints should all be negative in the optimum')
  %Nconstr*zopt - Kconstr';

  %transform back to original variables
  zopt = exp(zopt);
 
end
% now take the sensors out of this, and replace the values in x0 with these values
zopt = [fullxvector(zopt,BCs0)];
% check if we indeed have a critical point
disp('is dO/dx(z) == 0?')
dobjdx = dObjdx(zopt(:))
if norm(dobjdx) > 1e-2
  disp('not a local minimum. Exiting...');
  return;
end
