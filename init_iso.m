Globals1D;
GlobalsGR;

g00_exact = -power((2.*x-M)./(2.*x+M),2);
g01_exact = 0.*x;
g11_exact = power((2.*x+M)./(2.*x),4);
Pi00_exact = 0.*x;
Pi01_exact = 0.*x;
Pi11_exact = 0.*x;
Phi00_exact = 8.*(M-2.*x).*M./power(2.*x+M,3);
Phi01_exact = 0.*x;
Phi11_exact = -2.*M.*power((2.*x+M)./(2.*x),3)./power(x,2);
S_exact = 2.*log(1+M./(2.*x));
Pi_S_exact = 0.*x;
Phi_S_exact = -(M./power(x,2))./(1+M./(2.*x));
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

H0 = 0.*x;
H1 = 2./x./(1+M./(2.*x))./(1-M./(2.*x));

deriH00 = 0.*x;
deriH01 = 0.*x;
deriH10 = 0.*x;
deriH11 = -8.*(power(M,2)+4.*power(x,2))./power(power(M,2)-4.*power(x,2),2);
