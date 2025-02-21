function plotdata(t,x,y,mus,vs);

  figure(1)
  clf
  subplot(1,4,1)
  plot(t,x,'LineWidth',2)
  title('metab conc')

  subplot(1,4,2)
  plot(t,y,t,sum(y,2),'LineWidth',2)
  set(gca,'YLim',[0,1])
  legend('e1','e2','r1')
  title('enzymes and ribosomes')

  subplot(1,4,3)
  plot(t,mus,'LineWidth',2)
  legend('mu','Location','SouthEast')
  set(gca,'YLim',[0,1.2*max(mus(:))])
  title('growth rates') 

  subplot(1,4,4)
  plot(t,vs,'LineWidth',2)
  legend('fluxes','Location','SouthEast')
  set(gca,'YLim',[0,1.2*max(vs(:))])
  title('fluxes') 
 
