Globals1D;
GlobalsGR;

g00_exact = -(1-2./x);
g01_exact = sqrt(2./x);
g11_exact = 1+0.*x;
Pi00_exact = -2./power(x,2).*sqrt(2./x);
Pi01_exact = -1./power(x,2);
Pi11_exact = 0.*x;
Phi00_exact = -2./power(x,2);
Phi01_exact = -sqrt(1./(2.*x))./x;
Phi11_exact = 0.*x;
S_exact = 0.*x;
Pi_S_exact = 0.*x;
Phi_S_exact = 0.*x;
psi_exact = 0.*x;
Pi_psi_exact = 0.*x;
Phi_psi_exact = 0.*x;

g00 = g00_exact;
g01 = g01_exact;
g11 = g11_exact;
Pi00 = Pi00_exact;
Pi01 = Pi01_exact;
Pi11 = Pi11_exact;
Phi00 = Phi00_exact;
Phi01 = Phi01_exact;
Phi11 = Phi11_exact;
S = S_exact;
Pi_S = Pi_S_exact;
Phi_S = Phi_S_exact;
psi = psi_exact;
Pi_psi = Pi_psi_exact;
Phi_psi = Phi_psi_exact;

H0 = sqrt(1./(2.*x)).*(2+x)./power(x,2);
H1 = 1.*(1+0)./power(x,2) + 2./x;

deriH00 = 0.*x;
deriH01 = 0.*x;
deriH10 = -2.*sqrt(2).*(10+3.*x)./(8.*power(1./x, 1.5).*power(x, 5));
deriH11 = -2.*(0+1)./power(x,3) - 2./(x.*x);
