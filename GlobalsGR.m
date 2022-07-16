global g00 g01 g11 Pi00 Pi01 Pi11 Phi00 Phi01 Phi11 S Pi_S Phi_S psi Pi_psi Phi_psi
global rhs_g00 rhs_g01 rhs_g11 rhs_Pi00 rhs_Pi01 rhs_Pi11 rhs_Phi00 rhs_Phi01 rhs_Phi11 
global rhs_S rhs_Pi_S rhs_Phi_S rhs_psi rhs_Pi_psi rhs_Phi_psi
global lapse shift
global g00_exact g01_exact g11_exact Pi00_exact Pi01_exact Pi11_exact Phi00_exact Phi01_exact Phi11_exact
global S_exact Pi_S_exact Phi_S_exact psi_exact Pi_psi_exact Phi_psi_exact
global H0 H1 deriH00 deriH01 deriH10 deriH11
global DIRICHLET FREEZING bdry_type
%global T_scalar T00 T01 T11
%global invg00, invg01, invg11, lapse, shift, normal0, normal1, gamma11;
%global gamma000, gamma001, gamma011, gamma100, gamma101, gamma111, gamma0, gamma1;
global C0 C1 Cr00 Cr01 Cr11
global USE_S USE_psi
global AH_indicator;

M = 1.0;

%usual parameter
%paragamma0 = 3.*exp(-power(x./8, 2)./2)+0.1;
%paragamma1 = -1;
%paragamma2 = exp(-power(x./8, 2)./2)+0.1;
%paragamma4 = 0.5;
%paragamma5 = 0.5;

%testing paramater
paragamma0 = 1;
paragamma1 = -1;
paragamma2 = 1;
paragamma4 = 0;
paragamma5 = 0;

