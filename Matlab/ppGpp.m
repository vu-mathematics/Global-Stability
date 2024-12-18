function [a,p,m,r,mu,s] = ppGpp();

s=[0.01:0.1:1];
km=100; % kcat uptake
kr = 1; %kcat ribo
kR=75;kS=log(2)/30; % ppgpp kcats
Ks = 0.1; % Monod constant for substrate uptake
Ka = 100; % uptake inhibition by amino acid
Kr=.1; % aa saturation specificity
KR=.1; KS=.1; % affinity constants for ppGpp 
			%synth and degr by RelA and SpoT resp.
Kp=1; % sensitivity of ppGpp down regulation of ribo synth
Na=100;pT=1;


for i=1:length(s);
  pars = {s(i),km,kr,kR,kS,Ks,Ka,Kr,KR,KS,Kp,Na,pT};
  [t,y] = ppGppRun(10000,s(i),pars);
  y = y(end,:);
  [a(i),p(i),m(i),r(i)] = concentrations(y,pars);
  [f,v] = kinetics(a(i),p(i),m(i),r(i),pars);
  faa(i) = f(2)/Na;
  mu(i) = (v(5)+v(6))/(r(i)+m(i));
end

figure(1);
clf
subplot(2,2,1)
plot(mu,p,'-o','LineWidth',2);
xlabel('growth rate')
ylabel('ppGpp')

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



%%%%%%%%%%%%%%%%%%%%%%%

function [t,y] = ppGppRun(T,s,pars);



y0=[.01;0.2;1;1];

options = odeset('RelTol',1e-9,'AbsTol',1e-9);
PPGPP = @(t,y) ppgppODE(t,y,pars);
[t,y] = ode15s(PPGPP,[0,T],y0,options);




%%%%%%%%%%%%%%%%

function dy = ppgppODE(t,y,pars);

% extract the current metab, mRNA, protein concentrations from y
[a,p,m,r] = concentrations(y,pars);
% get kinetics and flux values.
[f,v] = kinetics(a,p,m,r,pars);

mu = (v(5)+v(6))/(r+m);

% compute dx/dt using a simple x' = Nv.
da = v(1)-v(2);
dp = v(3)-v(4);
dm = v(5)-mu*m;
dr = v(6)-mu*r;

dy = [da;dp;dm;dr];

%%%%%%%%%%

function [f,v] = kinetics(a,p,m,r,pars);

[s,km,kr,kR,kS,Ks,Ka,Kr,KR,KS,Kp,Na,pT] = pars{:};

fasynth = km * s/(s+Ks) * 1/(1+a/Ka);
facons = Na * kr * a/(a+Kr);

fpsynth = kR * Kr/(a+Kr) * (pT-p);%/(pT-p+KR);
fpdegr = kS * p;%/(p+KS);

frsynth = 1 / (1+p/Kp) * kr * a/(a+Kr);
fmsynth = (1- 1/ (1+p/Kp)) * kr * a/(a+Kr);


vasynth = m*fasynth;
vacons = r*facons;

vpsynth = fpsynth;
vpdegr = fpdegr;

vmsynth = r*fmsynth;
vrsynth = r*frsynth;

f = [fasynth;facons;fpsynth;fpdegr;fmsynth;frsynth];
v = [vasynth;vacons;vpsynth;vpdegr;vmsynth;vrsynth];

%%%%%%%%%%%%%

function [a,p,m,r] = concentrations(y,pars);

a = y(1);
p = y(2);
m = y(3);
r = y(4);
