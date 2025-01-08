function alphaopt=opt_egm(x0);

alpha = [0.5;0.3;0.2];


MU = @(alpha) -qorac(10000,alpha,x0,0);

A = [1 1 1;
     -1 -1 -1;	
     -1 0 0 ;	 
     0 -1 0 ; 
     0 0 -1];
B = [1;-1;0;0;0];

%opts = optimset('TolFun',1e-7);
opts = optimoptions(@fmincon,'Algorithm','sqp');
[alphaopt] = fmincon(MU,alpha,A,B,[],[],[],[],[],opts)
