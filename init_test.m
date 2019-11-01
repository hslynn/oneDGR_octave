Globals1D;
GlobalsGR;

g00_exact = -(1-2./power(x,2)-x);
g01_exact = 2./power(x,2)+x;
g11_exact = 1+3./power(x,2)+x;
Pi00_exact = -((8+power(x,3).*(2+3.*x+power(x,4)))./
     (power(x,3).*(3+power(x,2)+power(x,3)).*
       power((-2+power(x,2).*(1+x.*(-1+1.*x)))./
         (power(x,2).*(3+power(x,2)+power(x,3))),0.5)));
Pi01_exact = -((8+power(x,3).*(2+3.*x+power(x,4)))./
     (power(x,3).*(3+power(x,2)+power(x,3)).*
       power((-2+power(x,2).*(1+x.*(-1+1.*x)))./
         (power(x,2).*(3+power(x,2)+power(x,3))),0.5)));
Pi11_exact = -((12+power(x,3).*(4+3.*x+power(x,4)))./
     (power(x,3).*(3+power(x,2)+power(x,3)).*
       power((-2+power(x,2).*(1+x.*(-1+1.*x)))./
         (power(x,2).*(3+power(x,2)+power(x,3))),0.5)));
Phi00_exact = 1-4./power(x,3);
Phi01_exact = 1-4./power(x,3);
Phi11_exact = 1-6./power(x,3);

g00 = g00_exact;
g01 = g01_exact;
g11 = g11_exact;
Pi00 = Pi00_exact;
Pi01 = Pi01_exact;
Pi11 = Pi11_exact;
Phi00 = Phi00_exact;
Phi01 = Phi01_exact;
Phi11 = Phi11_exact;

H0 = 2./power(x,2);
H1 = 2.*(1+0)./power(x,2); %+ 2./x;

deriH00 = 0.*x;
deriH01 = 0.*x;
deriH10 = -4./power(x,3);
deriH11 = -2.*(0+2)./power(x,3); %- 2./(x.*x);
