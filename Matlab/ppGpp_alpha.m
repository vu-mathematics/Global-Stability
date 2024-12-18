function [mu,a,m,r,s] = ppGpp_alpha(s,alpha,doplot);

if nargin == 0
  s=[0.01:0.1:1];alpha=0.5;
elseif nargin == 1
  alpha=0.5;
end 
  
km = 1; % kcat uptake
kr = 1; %kcat ribo
Ks = 0.1; % Monod constant for substrate uptake
Ka = 1; % uptake inhibition by amino acid
Kr = 1; % aa saturation specificity
Na = 10; % number of aa's per protein


for i=1:length(s);
  pars = {s(i),km,kr,Ks,Ka,Kr,Na,alpha};
  [t,y] = ppGppRun(10000,s(i),pars);
  y = y(end,:);
  [a(i),m(i),r(i)] = concentrations(y,pars);
  [f,v] = kinetics(a(i),m(i),r(i),pars);
  faa(i) = f(2)/Na;
  mu(i) = (v(3)+v(4))/(r(i)+m(i));
end


if doplot
  figure(1);
  clf

  subplot(2,2,2)
  plot(s,mu,'o-','LineWidth',2);
  xlabel('nutrient conc')
  ylabel('growth rate')

  subplot(2,2,3)
  plot(mu,r./(r+m),'o-','LineWidth',2);
  xlabel('growth rate')
  ylabel('r/(m+r)')
 
  subplot(2,2,4)
  plot(mu,faa,'o-','LineWidth',2);
  xlabel('growth rate')
  ylabel('aa saturation')
end


%%%%%%%%%%%%%%%%%%%%%%%

function [t,y] = ppGppRun(T,s,pars);

y0=[.01;1;1];

options = odeset('RelTol',1e-9,'AbsTol',1e-9);
PPGPP = @(t,y) ppgppODE(t,y,pars);
[t,y] = ode15s(PPGPP,[0,T],y0,options);

%%%%%%%%%%%%%%%%

function dy = ppgppODE(t,y,pars);

% extract the current metab, mRNA, protein concentrations from y
[a,m,r] = concentrations(y,pars);
% get kinetics and flux values.
[f,v] = kinetics(a,m,r,pars);

mu = (v(3)+v(4))/(r+m);

% compute dx/dt using a simple x' = Nv.
da = v(1)-v(2);
dm = v(3)-mu*m;
dr = v(4)-mu*r;

dy = [da;dm;dr];

%%%%%%%%%%

function [f,v] = kinetics(a,m,r,pars);

[s,km,kr,Ks,Ka,Kr,Na,alpha] = pars{:};

fasynth = km * s/(s+Ks) * 1/(1+a/Ka);
facons = Na * kr * a/(a+Kr);

frsynth = alpha * kr * a/(a+Kr);
fmsynth = (1- alpha) * kr * a/(a+Kr);


vasynth = m*fasynth;
vacons = r*facons;

vmsynth = r*fmsynth;
vrsynth = r*frsynth;

f = [fasynth;facons;fmsynth;frsynth];
v = [vasynth;vacons;vmsynth;vrsynth];

%%%%%%%%%%%%%

function [a,m,r] = concentrations(y,pars);

a = y(1);
m = y(2);
r = y(3);
