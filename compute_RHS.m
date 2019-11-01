function [rhs_g00, rhs_g01, rhs_g11, rhs_Pi00, rhs_Pi01, rhs_Pi11, rhs_Phi00, rhs_Phi01, rhs_Phi11] = compute_RHS
Globals1D;
GlobalsGR;
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

gamma0 = invg00.*gamma000 + 2.*invg01.*gamma001 + invg11.*gamma011;
gamma1 = invg00.*gamma100 + 2.*invg01.*gamma101 + invg11.*gamma111;% - 2./x;

%constraints
C0 = H0+gamma0;
C1 = H1+gamma1;

%source terms
src_g00 = -lapse.*Pi00 - paragamma1.*shift.*Phi00;
src_g01 = -lapse.*Pi01 - paragamma1.*shift.*Phi01;
src_g11 = -lapse.*Pi11 - paragamma1.*shift.*Phi11;

src_Pi00 = (2.*lapse.*(invg00.*(gamma11.*Phi00.*Phi00-Pi00.*Pi00-invg00.*gamma000.*gamma000
        -2.*invg01.*gamma000.*gamma001 - invg11.*gamma001.*gamma001)
        + invg01.*(gamma11.*Phi00.*Phi01-Pi00.*Pi01-invg00.*gamma000.*gamma001
        - invg01.*(gamma000.*gamma011+gamma001.*gamma001)-invg11.*gamma001.*gamma011)
        + invg01.*(gamma11.*Phi00.*Phi01-Pi00.*Pi01-invg00.*gamma001.*gamma000
        - invg01.*(gamma001.*gamma001+gamma011.*gamma000)-invg11.*gamma011.*gamma001)
        + invg11.*(gamma11.*Phi01.*Phi01-Pi01.*Pi01-invg00.*gamma001.*gamma001
        - invg01.*2.*gamma001.*gamma011-invg11.*gamma011.*gamma011))
        %term 1
        - 0.5.*lapse.*Pi00.*(normal0.*normal0.*Pi00+2.*normal0.*normal1.*Pi01+normal1.*normal1.*Pi11)
        %term 2
        - lapse.*gamma11.*Phi00.*(normal0.*Pi01+normal1.*Pi11)
        %term 3
        - 2.*lapse.*(deriH00 + invg00.*gamma000.*(paragamma4.*C0-H0) + invg01.*gamma000.*(paragamma4.*C1-(H1+2./x))
        + invg01.*gamma100.*(paragamma4.*C0-H0) + invg11.*gamma100.*(paragamma4.*C1-(H1+2./x))
        - 0.5.*paragamma5.*g00.*(invg00.*gamma0.*C0+invg01.*(gamma0.*C1+gamma1.*C0)+invg11.*gamma1.*C1))
        %term 4
        + lapse.*paragamma0.*((-2.*lapse-g00.*normal0).*C0+(-g00.*normal1).*C1)
        %term 5
        - paragamma1.*paragamma2.*shift.*Phi00);
        %term 6
        %- 4.*lapse./power(x,2).*(gamma11.*g01.*-lapse.*(-normal1)).*(gamma11.*g01-lapse.*(-normal1)));
        %term 7


src_Pi01 = (2.*lapse.*(invg00.*(gamma11.*Phi00.*Phi01-Pi00.*Pi01-invg00.*gamma000.*gamma100 
        - invg01.*(gamma000.*gamma101+gamma001.*gamma100) - invg11.*gamma001.*gamma101)
        + invg01.*(gamma11.*Phi00.*Phi11-Pi00.*Pi11-invg00.*gamma000.*gamma101
        - invg01.*(gamma000.*gamma111+gamma001.*gamma101)-invg11.*gamma001.*gamma111)
        + invg01.*(gamma11.*Phi01.*Phi01-Pi01.*Pi01-invg00.*gamma001.*gamma100
        - invg01.*(gamma001.*gamma101+gamma011.*gamma100)-invg11.*gamma011.*gamma101)
        + invg11.*(gamma11.*Phi01.*Phi11-Pi01.*Pi11-invg00.*gamma001.*gamma101
        - invg01.*(gamma001.*gamma111+gamma011.*gamma101)-invg11.*gamma011.*gamma111))
        %term 1
        - 0.5.*lapse.*Pi01.*(normal0.*normal0.*Pi00+normal0.*normal1.*2.*Pi01+normal1.*normal1.*Pi11)
        %term 2
        - lapse.*gamma11.*Phi01.*(normal0.*Pi01+normal1.*Pi11) 
        %term 3
        - 2.*lapse.*(0.5.*(deriH01+deriH10)+invg00.*gamma001.*(paragamma4.*C0-H0)
        + invg01.*gamma001.*(paragamma4.*C1-(H1+2./x))
        + invg01.*gamma101.*(paragamma4.*C0-H0) + invg11.*gamma101.*(paragamma4.*C1-(H1+2./x))
        - 0.5.*paragamma5.*g01.*(invg00.*gamma0.*C0+invg01.*(gamma0.*C1+gamma1.*C0)+invg11.*gamma1.*C1))
        %term 4
        + lapse.*paragamma0.*((-g01.*normal0).*C0+(-lapse-g01.*normal1).*C1)
        %term 5
        - paragamma1.*paragamma2.*shift.*Phi01);
        %term 6
        %- 4.*lapse./power(x,2).*(gamma11.*g01-lapse.*(-normal1)).*(gamma11.*g11));
        %term 7


src_Pi11 = (2.*lapse.*(invg00.*(gamma11.*Phi01.*Phi01-Pi01.*Pi01-invg00.*gamma100.*gamma100
        -2.*invg01.*gamma100.*gamma101 - invg11.*gamma101.*gamma101)
        + invg01.*(gamma11.*Phi01.*Phi11-Pi01.*Pi11-invg00.*gamma100.*gamma101
        - invg01.*(gamma100.*gamma111+gamma101.*gamma101)-invg11.*gamma101.*gamma111)
        + invg01.*(gamma11.*Phi11.*Phi01-Pi11.*Pi01-invg00.*gamma101.*gamma100
        - invg01.*(gamma101.*gamma101+gamma111.*gamma100)-invg11.*gamma111.*gamma101)
        + invg11.*(gamma11.*Phi11.*Phi11-Pi11.*Pi11-invg00.*gamma101.*gamma101
        - invg01.*2.*gamma101.*gamma111-invg11.*gamma111.*gamma111))
        %term 1
        - 0.5.*lapse.*Pi11.*(normal0.*normal0.*Pi00+normal0.*normal1.*2.*Pi01+normal1.*normal1.*Pi11)
        %term 2
        - lapse.*gamma11.*Phi11.*(normal0.*Pi01+normal1.*Pi11)
        %term 3
        - 2.*lapse.*(deriH11 + invg00.*gamma011.*(paragamma4.*C0-H0) + invg01.*gamma011.*(paragamma4.*C1-(H1+2./x))
        + invg01.*gamma111.*(paragamma4.*C0-H0) + invg11.*gamma111.*(paragamma4.*C1-(H1+2./x))
        - 0.5.*paragamma5.*g11.*(invg00.*gamma0.*C0+invg01.*(gamma0.*C1+gamma1.*C0)+invg11.*gamma1.*C1))
        %term 4
        + lapse.*paragamma0.*((-g11.*normal0).*C0+(-g11.*normal1).*C1)
        %term 5
        - paragamma1.*paragamma2.*shift.*Phi11);
        %term 6
        %- 4.*lapse./power(x,2));
        %term 7

src_Phi00 = lapse.*(0.5.*Pi00.*(normal0.*normal0.*Phi00 + 2.*normal0.*normal1.*Phi01 + normal1.*normal1.*Phi11)
        + gamma11.*Phi00.*(normal0.*Phi01 + normal1.*Phi11) - paragamma2.*Phi00);
src_Phi01 = lapse.*(0.5.*Pi01.*(normal0.*normal0.*Phi00 + 2.*normal0.*normal1.*Phi01 + normal1.*normal1.*Phi11)
        + gamma11.*Phi01.*(normal0.*Phi01 + normal1.*Phi11) - paragamma2.*Phi01);
src_Phi11 = lapse.*(0.5.*Pi11.*(normal0.*normal0.*Phi00 + 2.*normal0.*normal1.*Phi01 + normal1.*normal1.*Phi11)
        + gamma11.*Phi11.*(normal0.*Phi01 + normal1.*Phi11) - paragamma2.*Phi11);

%Hhat terms
[deri_plus_g00, deri_minus_g00] = deri(g00, g00, g00);
[deri_plus_g01, deri_minus_g01] = deri(g01, g01, g01);
[deri_plus_g11, deri_minus_g11] = deri(g11, g11, g11);
[deri_plus_Pi00, deri_minus_Pi00] = deri(Pi00, Pi00, Pi00);
[deri_plus_Pi01, deri_minus_Pi01] = deri(Pi01, Pi01, Pi01);
[deri_plus_Pi11, deri_minus_Pi11] = deri(Pi11, Pi11, Pi11);
[deri_plus_Phi00, deri_minus_Phi00] = deri(Phi00, Phi00, Phi00);
[deri_plus_Phi01, deri_minus_Phi01] = deri(Phi01, Phi01, Phi01);
[deri_plus_Phi11, deri_minus_Phi11] = deri(Phi11, Phi11, Phi11);

max_alpha = 1.0;
avg_deri_g00 = 0.5.*(deri_plus_g00+deri_minus_g00);
avg_deri_g01 = 0.5.*(deri_plus_g01+deri_minus_g01);
avg_deri_g11 = 0.5.*(deri_plus_g11+deri_minus_g11);
avg_deri_Pi00 = 0.5.*(deri_plus_Pi00+deri_minus_Pi00);
avg_deri_Pi01 = 0.5.*(deri_plus_Pi01+deri_minus_Pi01);
avg_deri_Pi11 = 0.5.*(deri_plus_Pi11+deri_minus_Pi11);
avg_deri_Phi00 = 0.5.*(deri_plus_Phi00+deri_minus_Phi00);
avg_deri_Phi01 = 0.5.*(deri_plus_Phi01+deri_minus_Phi01);
avg_deri_Phi11 = 0.5.*(deri_plus_Phi11+deri_minus_Phi11);

Hhat_g00 = -(1+paragamma1).*shift.*avg_deri_g00 - max_alpha.*0.5.*(deri_plus_g00-deri_minus_g00);
Hhat_g01 = -(1+paragamma1).*shift.*avg_deri_g01 - max_alpha.*0.5.*(deri_plus_g01-deri_minus_g01);
Hhat_g11 = -(1+paragamma1).*shift.*avg_deri_g11 - max_alpha.*0.5.*(deri_plus_g11-deri_minus_g11);

Hhat_Pi00 = -paragamma1.*paragamma2.*shift.*avg_deri_g00 - shift.*avg_deri_Pi00 ...
        + lapse.*gamma11.*avg_deri_Phi00 - max_alpha.*0.5.*(deri_plus_Pi00-deri_minus_Pi00);
Hhat_Pi01 = -paragamma1.*paragamma2.*shift.*avg_deri_g01 - shift.*avg_deri_Pi01 ...
        + lapse.*gamma11.*avg_deri_Phi01 - max_alpha.*0.5.*(deri_plus_Pi01-deri_minus_Pi01);
Hhat_Pi11 = -paragamma1.*paragamma2.*shift.*avg_deri_g11 - shift.*avg_deri_Pi11 ...
        + lapse.*gamma11.*avg_deri_Phi11 - max_alpha.*0.5.*(deri_plus_Pi11-deri_minus_Pi11);

Hhat_Phi00 = -paragamma2.*lapse.*avg_deri_g00 + lapse.*avg_deri_Pi00 - shift.*avg_deri_Phi00 ...
        - max_alpha.*0.5.*(deri_plus_Phi00-deri_minus_Phi00);
Hhat_Phi01 = -paragamma2.*lapse.*avg_deri_g01 + lapse.*avg_deri_Pi01 - shift.*avg_deri_Phi01 ...
        - max_alpha.*0.5.*(deri_plus_Phi01-deri_minus_Phi01);
Hhat_Phi11 = -paragamma2.*lapse.*avg_deri_g11 + lapse.*avg_deri_Pi11 - shift.*avg_deri_Phi11 ...
        - max_alpha.*0.5.*(deri_plus_Phi11-deri_minus_Phi11);

%RHS
rhs_g00 = src_g00 - Hhat_g00;
rhs_g01 = src_g01 - Hhat_g01;
rhs_g11 = src_g11 - Hhat_g11;
rhs_Pi00 = src_Pi00 - Hhat_Pi00;
rhs_Pi01 = src_Pi01 - Hhat_Pi01;
rhs_Pi11 = src_Pi11 - Hhat_Pi11;
rhs_Phi00 = src_Phi00 - Hhat_Phi00;
rhs_Phi01 = src_Phi01 - Hhat_Phi01;
rhs_Phi11 = src_Phi11 - Hhat_Phi11;

%constraints
Cr11 = avg_deri_g11 - Phi11;


return
