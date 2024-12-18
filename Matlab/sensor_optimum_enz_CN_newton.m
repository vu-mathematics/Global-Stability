function [zopt,Eopt,Jopt,Esatopt] = sensor_optimum_enz_CN;

% function [zopt,Eopt,Jopt,Esatopt] = sensor_optimum_enz_CN;
%
% Belongs to "Dynamic and robust adaptive control of optimal
% metabolic flux" by Planque et al. 
% 
% computes the optimal input-output relations for the metabolic pathway 
% in Figures 2 and 3.
%
% The script computes optima for varying external concentrations, by
% solving the dO/dx=0 equations using Newton's method, and then plots the enzyme synthesis rates
% against the sensor values in the optima.

% the following file contains the kinetics functions and their Jacobians.
% To make this file, first run the daes_CN.m file. It produces
% the CN_1.mat file, which can then be loaded here.
load('CN_1.mat');

params{1} = A;
params{3} = n;
params{5} = sensors;
params{7} = nenz;
params{8} = OUT;
params{9} = NN;

% we need a start optimum to continue 
BCs0=5;
zopt = [];
zopt(:,1) = compute_optimum(BCs0,params,Objsym,ftns,dObjdx,dObjdxx_newton);
Eopt(:,1) = Obj_parts(zopt(:,1)); 
som=sum(Eopt(:,1));
Eopt(:,1) = Eopt(:,1)/som;
Zopt = zopt;
zoptlast = zopt(:,1);

% compute from BC=5 to BC=1, so going down...
No=50;
BCs=linspace(5,.4,No);
i=2; notnan=1;
while (i<=No && notnan == 1)
  BC = BCs(i)
  zopt(:,i) = compute_optimum(BC,params,Objsym,ftns,dObjdx,dObjdxx_newton,zoptlast(2:4));
  if isnan(zopt(:,i))
    zopt = zopt(:,1:end-1);
    notnan=0;
  else
    zoptlast = zopt(:,i);
    Eopt(:,i) = Obj_parts(zopt(:,i)); som=sum(Eopt(:,i));
    Eopt(:,i) = Eopt(:,i)/som;
  end
  i = i+1;
end
zoptlow = zopt;
Eoptlow = Eopt;
BCslow=BCs(1:size(zopt,2));

% now start at BC=5 again, but going up
zopt = Zopt;
Eopt = [];
BCs=linspace(5,50,No);
i=2; notnan=1;
while (i<=No && notnan == 1) % we compute until the equations give NaN as result.
  BC = BCs(i);
  zopt(:,i) = compute_optimum(BC,params,Objsym,ftns,dObjdx,dObjdxx_newton,zopt(2:4,i-1));
  if isnan(zopt(:,i))
    zopt = zopt(:,1:end-1);
    notnan=0;
  else
    % compute enzyme synthesis functions == input-output relations
    Eopt(:,i) = Obj_parts(zopt(:,i)); 	
    som=sum(Eopt(:,i));
    Eopt(:,i) = Eopt(:,i)/som;
  end
  i = i+1;
end
zopthigh = zopt;
BCshigh=BCs(1:size(zopt,2));
Eopthigh = Eopt;

% collect everything together
zopt = [zoptlow zopthigh(:,2:end)];
BCs = [BCslow BCshigh(:,2:end)];
Eopt = [Eoptlow Eopthigh(:,2:end)];
[BCs,v] = sort(BCs);
zopt = zopt(:,v);
Eopt = Eopt(:,v);

% now compute the input-output relations
for i=1:length(Eopt)
  Fopt = ftns(zopt(:,i));
  Jopt(i) = Eopt(end,i)*Fopt(end)';
  Esatopt(i) = Fopt(end);
end 

% and plot the result
figure(1)
clf

%subplot(1,3,2)
hold on
p1 = plot(zopt(2,:),Eopt(1,:),'-','LineWidth',2);
p2 = plot(zopt(2,:),Eopt(2,:),'-','LineWidth',2);
p3 = plot(zopt(2,:),Eopt(3,:),'-','LineWidth',2);
p4 = plot(zopt(2,:),Eopt(4,:),'-','LineWidth',2);
hold off
legend([p1(1) p2(1) p3(1) p4(1)],'Enzyme 1', 'Enzyme 2', 'Enzyme 3', 'Enzyme 4')

print_fig(1,'pics/CN_1_input_output.jpg',[0 0 6 4]);


%%%%%%%%%%%%%%%%%%%%%%%%

function zopt = compute_optimum(BCs0,params,Objsym,ftns,dObjdx,dObjdxx_newton,z0);

% computes an optimum for given external concentrations using Newton's method on the dO/dx=0 equations. 
% See daes_CN.m for more details. 

A = params{1};
n = params{3};
sensors = params{5};
nenz = params{7};
OUT = params{8};
NN = params{9};

fullxvector=@(x,BCs) [BCs(1) x(1:end)'];

f = @(x) dObjdx(fullxvector(x,BCs0))';
df = @(x) dObjdxx_newton(fullxvector(x,BCs0));

if nargin < 7
    z0 = [    1.6430    1.1094    0.1865]'; % just some random IC that is in the admissable set.
end
zopt = fsolve(f,z0);

zopt = [fullxvector(zopt,BCs0)];
dobjdx = dObjdx(zopt(:));
if norm(dobjdx) > 1e-2
  zopt = nan;
  disp('not a local minimum. Exiting...');
  return;
end

