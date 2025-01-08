function [mu,t,xs,y,mus,vs] = qorac(T,alpha,x0s,doplot)

% function [mu,t,y,mus] = qorac(T,alpha,x0s,doplot) 
%
% the qORAc script for qORAC control of the simple metabolic
% network shown in Figure 1 of 
% 
% Planqu\'e et al, Sensing cellular growth rate 
% facilitates robust optimal adaptation to changing conditions
% 
% See the Supplementary Information, section S3 for additional
% explanation of these simulations.
%
% Parameters:
% T: time span of one run (t in [0,T])
% alpha: ribosome allocation fractions 
%	alpha = 0 means alpha is not supplied. a default is chosen further down
%   alpha = 1 means we do qORAC control
%   alpha = a vector means we chose a fixed allocation. Used in finding optimal 
%			allocations by varying alpha
% x0s: one or more nutrient concentrations. When x0s is a vector, the simulation
% 	does a qORAC run from one nutrient concentration to the next, with the 
% 	network adapting to these changes
% doplot: flag for plotting (doplot=1, useful for individual runs) or not (doplot=0, useful 
% 			when finding optimal allocations, so this script is called many times in one run
%			run of scan_x0)
% fixedalpha is a flag with which we can change from finding optima 
% 	(in which we need to find the growth rate for fixed alphas all the time)
% 	to one in which we do qORAC control, changing the nutrient concentrations
% 	once every so often.
%   this flag is a separate one from the alpha parameter above, since alpha is sometimes 
% 	overwritten (it must become a vector at some point in the simulation).


% the optimal allocations (alphas, also termed epsilon in the paper, 
% since they differ up to a scalar constant)
% are computed using scan_x0.m
load('alphaopts1.mat')
load('muopts1.mat')
load('x0opts');

% the starting enzyme concentrations are not optimal. take something random away from optimal
% we use the alphas for that, but the alphas are the same as the enzyme concentrations,
% up to a constant. 
%alpha0 = interp1(x0opts',alphaopt',x0s(1)) + 0.1*rand(1,3);
%if isnan(alpha0)
%  alpha0 = rand(1,3);
%end
%alpha0 = alpha0/(sum(alpha0));	% renormalise what you chose to sum = 1. 

% all parameters that are passed around in all functions below. 
% kcat1 and kcat2, kcat_ribo, current nutrient concentration x0, alpha, the optimal alphas and mu's, 
% and the qssa flag
pars = {[5,5], 0.3, x0s(1), alpha,alphaopt,muopts}; 

% just some initial condition
y0 = [0.2 0.3 0.5];
% set up tolerances for the ODE solver, that are stricter than the built in defaults.
options = odeset('RelTol',1e-8,'AbsTol',1e-8);
% define the right hand side of the enzyme dynamics equations
RHS = @(t,y) rhs(t,y,pars);
% then run the simulation
[t,y] = ode15s(RHS,[0,T],y0,options);
X0s = x0s(1)*ones(length(t),1);

% the remaining part of the script is to compute statistics for all the 
% solutions, such growth rates   
mus = zeros(length(t),1);
vs = zeros(length(t),3);
xs = zeros(length(t),2);



% take away the semicolons if you wish to inspect these variables during output
alpha(:,:);
mu = mus(end);


if length(x0s)>1
  mumax = interp1(x0opts',muopts,x0s(1));
  MUmaxs = mumax*ones(length(t),1);


  % for each nutrient concentration, we rerun the simulation above,
  % using the end of the previous run as the start of the next one
  for i=2:length(x0s)
    %update the nutrient concentration
    pars{3}=x0s(i);
    
	% (note that x0 has changed inside pars)
    RHS = @(t,y) rhs(t,y,pars);

    % take the end of the last run as the start of the new one
    y0 = y(end,:);
    Tend = t(end);
    % run the next bit, and add it to the previous one
    [t2,y2] = ode15s(RHS,[Tend,Tend+T],y0,options);

    t = [t;t2];
    y = [y;y2];
    X0s = [X0s; x0s(i)*ones(length(t2),1)];
	mumax = interp1(x0opts',muopts,x0s(i));
	MUmaxs = [MUmaxs; mumax*ones(length(t2),1)];

  end
else
  MUmaxs = mu*ones(length(t),1);
end

for i=1:length(t)
  X0 = [0.1; 0.2];
  pars{7} = y(i,:);
  pars{3} = X0s(i);
  QSSA = @(tt,x) qssa(tt,x,pars);
  [tx,x] = ode15s(QSSA,[0,1000],X0,options);
  xs(i,:) = x(end,1:2);
  [vs(i,:),~,mus(i)] = metab(xs(i,:),y(i,:),pars);
end
  
  
% plot routines. turn off plotting by setting doplot=0  
if doplot
  plotdata(t,xs,y,mus,vs);
else
  [alpha(:)'  mu]
end

% save the results to disk
data = [t xs y mus X0s MUmaxs];
save('qorac_example.txt','-ASCII','data');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dx = qssa(t,x,pars);

% a version of the rhs script used just for the qssa

% pull out all the parameters
[kr,kc,x0,alpha,alphaopts1,muopts1,enz] = pars{:};

% get the current metabolite network state (fluxes, metabolite concentrations, growth rate)
[v,dx,mu] = metab(x,enz,pars);

dx = dx(:);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function dy = rhs(t,y,pars);

% the main script in which the qORAC controlled network is defined

% pull out all the parameters
[kr,kc,x0,alpha,alphaopts1,muopts1] = pars{:};

pars{7} = y;
X0 = [0.1; 0.2];
options = odeset('RelTol',1e-8,'AbsTol',1e-8);
QSSA = @(t,x) qssa(t,x,pars);
[tx,x] = ode15s(QSSA,[0,1000],X0,options);
x = x(end,1:2);

% get the current metabolite network state (fluxes, metabolite concentrations, growth rate)
[v,~,mu] = metab(x,y,pars);

% if alpha is 1, do the qORAC control  
if alpha==1
  % as explained in the Supplementary Information, all optima are
  % precomputed, and we just pull out the one we want using linear interpolation on mu
  alpha = interp1(muopts1',alphaopts1',mu)
end  

e1 = y(1);  e2 = y(2);  r = y(3);

% the enzyme equations are simply synthesis minus dilution. 
de1 = v(3) * alpha(1) - mu * e1;
de2 = v(3) * alpha(2) - mu * e2;
dr  = v(3) * alpha(3) - mu * r;

% collect all d/dt's into one vector
dy = [de1;de2;dr];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [v,dx,mu] = metab(x,y,pars);

% script to obtain the state of the metabolic network 
% we actually do not do the QSSA in each time step; the metabolite dynamics are
% much faster than the enzyme changes, anyway, so the difference between
% enforcing QSSA and not are minimal.

% pull out all the parameters
[kcs,kr,x0,alpha,alphaopts1,muopts1] = pars{:};

% the kcats of the enzymatic reaction
kc1=kcs(1);
kc2=kcs(2);

% the metabolite concentrations
x1 = x(1); x2 = x(2);
% the enzyme / ribosome concentrations
e1 = y(1); e2 = y(2); r = y(3);

% the reaction fluxes
v(1) = kc1 * e1 * (x0-0.1*x1) /(1+2*x0+1*x1);
v(2) = kc2 * e2 * (x1-0.1*x2) / (1+x2+1*x1);
v(3) = kr * r * x2 / (.1+5*x2);

% total enzyme concentration
et = e1 + e2 + r;

% the growth rate is protein synthesis rate per unit protein
mu = v(3)  / et;

% the actual dynamical system for metabolism, with stoichiometry integrated.
dx(1) = v(1) - v(2) ;
dx(2) = v(2) - 10*v(3);

