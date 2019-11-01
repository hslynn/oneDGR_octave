function [rhs_Pi00_L2, rhs_Pi01_L2, rhs_Pi11_L2] = compute_RHS_stiff(g00, g01, g11, Phi00, Phi01, Phi11)
Globals1D;
%inverse metric
invg00 = g11./(g00.*g11-g01.*g01);
invg01 = -g01./(g00.*g11-g01.*g01);
invg11 = g00./(g00.*g11-g01.*g01);

%some auxiliary variabls
lapse = 1.0./power(-invg00, 0.5);
shift = -invg01./invg00;
normal0 = 1.0./lapse;
normal1 = -shift./lapse;
gamma11 =  1.0./g11;

%connections
gamma000 = gamma11.*0.5.*(2.*g01.*Phi00) - 0.5.*gamma11.*g01.*Phi00 - lapse.*Pi00 + 0.5.*lapse.*Pi00;
gamma001 = gamma11.*0.5.*(g01.*Phi01+g11.*Phi00) - gamma11.*0.5.*g01.*Phi01;
gamma011 = gamma11.*0.5.*(2.*g11.*Phi01) - 0.5.*gamma11.*g01.*Phi11 + 0.5.*lapse.*Pi11;

gamma100 = gamma11.*0.5.*(2.*g01.*Phi01) - 0.5.*Phi00 - lapse.*Pi01;
gamma101 = gamma11.*0.5.*(g01.*Phi11 + g11.*Phi01) - 0.5.*Phi01 - 0.5.*lapse.*Pi11;
gamma111 = gamma11.*0.5.*(2.*g11.*Phi11) - 0.5.*Phi11;

%source terms
src_g00 = -lapse.*Pi00 - paragamma1.*shift.*Phi00;
src_g01 = -lapse.*Pi01 - paragamma1.*shift.*Phi01;
src_g11 = -lapse.*Pi11 - paragamma1.*shift.*Phi11;

src_Pi00_L2 = - 2.*lapse.*(deriH00 + invg00.*gamma000.*(paragamma4.*C0-H0) + invg01.*gamma000.*(paragamma4.*C1-(H1+2./x))
        + invg01.*gamma100.*(paragamma4.*C0-H0) + invg11.*gamma100.*(paragamma4.*C1-(H1+2./x))
        - 0.5.*paragamma5.*g00.*(invg00.*gamma0.*C0+invg01.*(gamma0.*C1+gamma1.*C0)+invg11.*gamma1.*C1));

src_Pi01_L2 = - 2.*lapse.*(0.5.*(deriH01+deriH10)+invg00.*gamma001.*(paragamma4.*C0-H0)
        + invg01.*gamma001.*(paragamma4.*C1-(H1+2./x))
        + invg01.*gamma101.*(paragamma4.*C0-H0) + invg11.*gamma101.*(paragamma4.*C1-(H1+2./x))
        - 0.5.*paragamma5.*g01.*(invg00.*gamma0.*C0+invg01.*(gamma0.*C1+gamma1.*C0)+invg11.*gamma1.*C1));

src_Pi11_L2 = - 2.*lapse.*(deriH11 + invg00.*gamma011.*(paragamma4.*C0-H0) + invg01.*gamma011.*(paragamma4.*C1-(H1+2./x))
        + invg01.*gamma111.*(paragamma4.*C0-H0) + invg11.*gamma111.*(paragamma4.*C1-(H1+2./x))
        - 0.5.*paragamma5.*g11.*(invg00.*gamma0.*C0+invg01.*(gamma0.*C1+gamma1.*C0)+invg11.*gamma1.*C1));

rhs_Pi00_L2 = src_Pi00_L2;
rhs_Pi01_L2 = src_Pi01_L2;
rhs_Pi11_L2 = src_Pi11_L2;

return
