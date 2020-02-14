Globals1D;
GlobalsGR;

g00_exact = -(1-2.*M./x);
g01_exact = 2.*M./x;
g11_exact = 1+2.*M./x;
Pi00_exact = -2.*M./(x+2.*M).*sqrt(1+2.*M./x).*2.*M./power(x,2);
Pi01_exact = -2.*M./(x+2.*M).*sqrt(1+2.*M./x).*2.*M./power(x,2);
Pi11_exact = -2.*M./(x+2.*M).*sqrt(1+2.*M./x).*2.*M./power(x,2);
Phi00_exact = -2.*M./power(x,2);
Phi01_exact = -2.*M./power(x,2);
Phi11_exact = -2.*M./power(x,2);
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

H0 = 2.*M./power(x,2);
H1 = 2.*(M+0)./power(x,2) + 2./x;

deriH00 = 0.*x;
deriH01 = 0.*x;
deriH10 = -4.*M./power(x,3);
deriH11 = -2.*(0+2.*M)./power(x,3) - 2./(x.*x);
