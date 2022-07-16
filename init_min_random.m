Globals1D;
GlobalsGR;

g00_exact = -1+0.*x;
g01_exact = 0.*x;
g11_exact = 1+0.*x;
Pi00_exact = 0.*x;
Pi01_exact = 0.*x;
Pi11_exact = 0.*x;
Phi00_exact = 0.*x;
Phi01_exact = 0.*x;
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

H0 = 0.*x;
H1 = 0.*x + 2./x;

deriH00 = 0.*x;
deriH01 = 0.*x;
deriH10 = 0.*x;
deriH11 = 0.*x - 2./(x.*x);

%generate random perturbation
epsilon = power(10, -6);
g00 = g00 + (2.*rand(size(x)(1), size(x)(2))-1).*epsilon;
g01 = g01 + (2.*rand(size(x)(1), size(x)(2))-1).*epsilon;
g11 = g11 + (2.*rand(size(x)(1), size(x)(2))-1).*epsilon;
Pi00 = Pi00 + (2.*rand(size(x)(1), size(x)(2))-1).*epsilon;
Pi01 = Pi01 + (2.*rand(size(x)(1), size(x)(2))-1).*epsilon;
Pi11 = Pi11 + (2.*rand(size(x)(1), size(x)(2))-1).*epsilon;
Phi00 = Phi00 + (2.*rand(size(x)(1), size(x)(2))-1).*epsilon;
Phi01 = Phi01 + (2.*rand(size(x)(1), size(x)(2))-1).*epsilon;
Phi11 = Phi11 + (2.*rand(size(x)(1), size(x)(2))-1).*epsilon;
