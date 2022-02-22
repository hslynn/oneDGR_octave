function compute_RHS
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

gamma0 = invg00.*gamma000 + 2.*invg01.*gamma001 + invg11.*gamma011 ...
        - 2./x.*(gamma11.*g01.*(x.*Phi_S+1)-lapse.*(x.*Pi_S-normal1));
gamma1 = invg00.*gamma100 + 2.*invg01.*gamma101 + invg11.*gamma111 - 2./x.*(x.*Phi_S+1);

%energy momentum tensor
T_scalar = Pi_psi.*Pi_psi - gamma11.*Phi_psi.*Phi_psi;
T00 = gamma11.*gamma11.*g01.*g01.*Phi_psi.*Phi_psi + 2.*gamma11.*g01.*(-lapse).*Phi_psi.*Pi_psi ...
    + lapse.*lapse.*Pi_psi.*Pi_psi + 0.5.*g00.*T_scalar;
T01 = gamma11.*gamma11.*g01.*g11.*Phi_psi.*Phi_psi + gamma11.*g11.*(-lapse).*Phi_psi.*Pi_psi ...
    + 0.5.*g01.*T_scalar;
T11 = gamma11.*gamma11.*g11.*g11.*Phi_psi.*Phi_psi + 0.5.*g11.*T_scalar;


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
        - 2.*lapse.*(deriH00 + invg00.*gamma000.*(paragamma4.*C0-H0) + invg01.*gamma000.*(paragamma4.*C1-H1)
        + invg01.*gamma100.*(paragamma4.*C0-H0) + invg11.*gamma100.*(paragamma4.*C1-H1)
        - 0.5.*paragamma5.*g00.*(invg00.*gamma0.*C0+invg01.*(gamma0.*C1+gamma1.*C0)+invg11.*gamma1.*C1))
        %term 4
        + lapse.*paragamma0.*((-2.*lapse-g00.*normal0).*C0+(-g00.*normal1).*C1)
        %term 5
        - paragamma1.*paragamma2.*shift.*Phi00
        %term 6
        - 4.*lapse./power(x,2).*(gamma11.*g01.*(x.*Phi_S+1)-lapse.*(x.*Pi_S-normal1)).*(gamma11.*g01.*(x.*Phi_S+1)
        -lapse.*(x.*Pi_S-normal1)) 
        %term 7
        + 16.*pi.*lapse.*(T00-0.5.*T_scalar.*g00));

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
        + invg01.*gamma001.*(paragamma4.*C1-H1)
        + invg01.*gamma101.*(paragamma4.*C0-H0) + invg11.*gamma101.*(paragamma4.*C1-H1)
        - 0.5.*paragamma5.*g01.*(invg00.*gamma0.*C0+invg01.*(gamma0.*C1+gamma1.*C0)+invg11.*gamma1.*C1))
        %term 4
        + lapse.*paragamma0.*((-g01.*normal0).*C0+(-lapse-g01.*normal1).*C1)
        %term 5
        - paragamma1.*paragamma2.*shift.*Phi01
        %term 6
        - 4.*lapse./power(x,2).*(gamma11.*g01.*(x.*Phi_S+1)-lapse.*(x.*Pi_S-normal1)).*(gamma11.*g11.*(x.*Phi_S+1))
        %term 7
        + 16.*pi.*lapse.*(T01-0.5.*T_scalar.*g01));


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
        - 2.*lapse.*(deriH11 + invg00.*gamma011.*(paragamma4.*C0-H0) + invg01.*gamma011.*(paragamma4.*C1-H1)
        + invg01.*gamma111.*(paragamma4.*C0-H0) + invg11.*gamma111.*(paragamma4.*C1-H1)
        - 0.5.*paragamma5.*g11.*(invg00.*gamma0.*C0+invg01.*(gamma0.*C1+gamma1.*C0)+invg11.*gamma1.*C1))
        %term 4
        + lapse.*paragamma0.*((-g11.*normal0).*C0+(-g11.*normal1).*C1)
        %term 5
        - paragamma1.*paragamma2.*shift.*Phi11
        %term 6
        - 4.*lapse./power(x,2).*(gamma11.*g11.*(x.*Phi_S+1)).*(gamma11.*g11.*(x.*Phi_S+1))
        %term 7
        + 16.*pi.*lapse.*(T11-0.5.*T_scalar.*g11));


src_Phi00 = lapse.*(0.5.*Pi00.*(normal0.*normal0.*Phi00 + 2.*normal0.*normal1.*Phi01 + normal1.*normal1.*Phi11)
        + gamma11.*Phi00.*(normal0.*Phi01 + normal1.*Phi11) - paragamma2.*Phi00);
src_Phi01 = lapse.*(0.5.*Pi01.*(normal0.*normal0.*Phi00 + 2.*normal0.*normal1.*Phi01 + normal1.*normal1.*Phi11)
        + gamma11.*Phi01.*(normal0.*Phi01 + normal1.*Phi11) - paragamma2.*Phi01);
src_Phi11 = lapse.*(0.5.*Pi11.*(normal0.*normal0.*Phi00 + 2.*normal0.*normal1.*Phi01 + normal1.*normal1.*Phi11)
        + gamma11.*Phi11.*(normal0.*Phi01 + normal1.*Phi11) - paragamma2.*Phi11);

if (USE_S)
    src_S = -lapse.*Pi_S - paragamma1.*shift.*Phi_S;
    src_Pi_S = (-lapse.*Phi_S.*gamma11.*(normal0.*Pi01 + normal1.*Pi11) - 0.5.*lapse.*Pi_S.*(normal0.*normal0.*Pi00 
            + 2.*normal0.*normal1.*Pi01 + normal1.*normal1.*Pi11) - paragamma1.*paragamma2.*shift.*Phi_S 
            %row 1
            - 2.*lapse./power(x,2).*invg00.*(gamma11.*g01.*(x.*Phi_S+1)-lapse.*(x.*Pi_S-normal1)).*(gamma11.*g01.*(x.*Phi_S+1) 
            - lapse.*(x.*Pi_S-normal1)) 
            - 4.*lapse./power(x,2).*invg01.*(gamma11.*g01.*(x.*Phi_S+1)-lapse.*(x.*Pi_S-normal1)).*gamma11.*g11.*(x.*Phi_S+1) 
            - 2.*lapse./power(x,2).*invg11.*(gamma11.*g11.*(x.*Phi_S+1)).*(gamma11.*g11.*(x.*Phi_S+1)) 
            %row 2
            - lapse./x.*invg00.*H0.*(gamma11.*g01.*(x.*Phi_S+1)-lapse.*(x.*Pi_S-normal1)) 
            - lapse./x.*invg01.*H0.*gamma11.*g11.*(x.*Phi_S+1)-lapse./x.*invg01.*H1.*(gamma11.*g01.*(x.*Phi_S+1) 
            - lapse.*(x.*Pi_S-normal1)) 
            - lapse./x.*invg11.*H1.*gamma11.*g11.*(x.*Phi_S+1) 
            %row 3
            - 2.*lapse.*(Pi_S.*Pi_S-gamma11.*Phi_S.*Phi_S) + 4.*lapse./x.*(normal1.*Pi_S+gamma11.*Phi_S) 
            + 3./power(x,2).*(lapse.*gamma11+shift.*normal1) + lapse./power(x,2).*exp(-2.*S)
            %row 4
            %+ 8./x.*lapse.*Phi_S
            - 0.5.*lapse.*paragamma0.*(normal0.*C0+normal1.*C1));
            %TBD
    src_Phi_S = lapse.*(0.5.*Pi_S.*(normal0.*normal0.*Phi00 + 2.*normal0.*normal1.*Phi01 + normal1.*normal1.*Phi11) 
            + gamma11.*Pi_S.*(normal0.*Phi01 + normal1.*Phi11) - paragamma2.*Phi_S);
end;

if (USE_psi)
    src_psi = -lapse.*Pi_psi - paragamma1.*shift.*Phi_psi;
    src_Pi_psi = lapse.*gamma11.*gamma1.*Phi_psi + lapse.*Pi_psi.*(invg00.*gamma0.*(-lapse)+invg01.*gamma1.*(-lapse)) ...
        - paragamma1.*paragamma2.*shift.*Phi_psi - lapse.*(gamma11.*Phi_psi.*(normal0.*Pi01+normal1.*Pi11) ...
        + 0.5.*Pi_psi.*(normal0.*normal0.*Pi00 + 2.*normal0.*normal1.*Pi01 + normal1.*normal1.*Pi11));
    src_Phi_psi = lapse.*(0.5.*Pi_psi.*(normal0.*normal0.*Phi00 + 2.*normal0.*normal1.*Phi01 + normal1.*normal1.*Phi11) ...
        + gamma11.*Phi_psi.*(normal0.*Phi01 + normal1.*Phi11) - paragamma2.*Phi_psi);
end;



%Hhat terms
if (bdry_type == DIRICHLET)
    g00_rbdry = g00_exact(vmapO);
    g01_rbdry = g01_exact(vmapO);
    g11_rbdry = g11_exact(vmapO);
    Pi00_rbdry = Pi00_exact(vmapO);
    Pi01_rbdry = Pi01_exact(vmapO);
    Pi11_rbdry = Pi11_exact(vmapO);
    Phi00_rbdry = Phi00_exact(vmapO);
    Phi01_rbdry = Phi01_exact(vmapO);
    Phi11_rbdry = Phi11_exact(vmapO);
    if (USE_S)
        S_rbdry = S_exact(vmapO);
        Pi_S_rbdry = Pi_S_exact(vmapO);
        Phi_S_rbdry = Phi_S_exact(vmapO);
    end;
    if (USE_psi)
        psi_rbdry = psi_exact(vmapO);
        Pi_psi_rbdry = Pi_psi_exact(vmapO);
        Phi_psi_rbdry = Phi_psi_exact(vmapO);
    end;
else
    g00_rbdry = g00(vmapO);
    g01_rbdry = g01(vmapO);
    g11_rbdry = g11(vmapO);
    Pi00_rbdry = Pi00(vmapO);
    Pi01_rbdry = Pi01(vmapO);
    Pi11_rbdry = Pi11(vmapO);
    Phi00_rbdry = Phi00(vmapO);
    Phi01_rbdry = Phi01(vmapO);
    Phi11_rbdry = Phi11(vmapO);
    if (USE_S)
        S_rbdry = S(vmapO);
        Pi_S_rbdry = Pi_S(vmapO);
        Phi_S_rbdry = Phi_S(vmapO);
    end;    
    if (USE_psi)
        psi_rbdry = psi(vmapO);
        Pi_psi_rbdry = Pi_psi(vmapO);
        Phi_psi_rbdry = Phi_psi(vmapO);
    end;
end;

if (inB < 2)
    g00_lbdry = g00(vmapI);
    g01_lbdry = g01(vmapI);
    g11_lbdry = g11(vmapI);
    Pi00_lbdry = Pi00(vmapI);
    Pi01_lbdry = Pi01(vmapI);
    Pi11_lbdry = Pi11(vmapI);
    Phi00_lbdry = Phi00(vmapI);
    Phi01_lbdry = Phi01(vmapI);
    Phi11_lbdry = Phi11(vmapI);
    if (USE_S)
        S_lbdry = S(vmapI);
        Pi_S_lbdry = Pi_S(vmapI);
        Phi_S_lbdry = Phi_S(vmapI);
    end;
    if (USE_psi)
        psi_lbdry = psi(vmapI);
        Pi_psi_lbdry = Pi_psi(vmapI);
        Phi_psi_lbdry = Phi_psi(vmapI);
    end;
else
    g00_lbdry = g00_exact(vmapI);
    g01_lbdry = g01_exact(vmapI);
    g11_lbdry = g11_exact(vmapI);
    Pi00_lbdry = Pi00_exact(vmapI);
    Pi01_lbdry = Pi01_exact(vmapI);
    Pi11_lbdry = Pi11_exact(vmapI);
    Phi00_lbdry = Phi00_exact(vmapI);
    Phi01_lbdry = Phi01_exact(vmapI);
    Phi11_lbdry = Phi11_exact(vmapI);
    if (USE_S)
        S_lbdry = S_exact(vmapI);
        Pi_S_lbdry = Pi_S_exact(vmapI);
        Phi_S_lbdry = Phi_S_exact(vmapI);
    end;
    if (USE_psi)
        psi_lbdry = psi_exact(vmapI);
        Pi_psi_lbdry = Pi_psi_exact(vmapI);
        Phi_psi_lbdry = Phi_psi_exact(vmapI);
    end;
end;


[deri_plus_g00, deri_minus_g00] = deri(g00, g00_lbdry, g00_rbdry);
[deri_plus_g01, deri_minus_g01] = deri(g01, g01_lbdry, g01_rbdry);
[deri_plus_g11, deri_minus_g11] = deri(g11, g11_lbdry, g11_rbdry);
[deri_plus_Pi00, deri_minus_Pi00] = deri(Pi00, Pi00_lbdry, Pi00_rbdry);
[deri_plus_Pi01, deri_minus_Pi01] = deri(Pi01, Pi01_lbdry, Pi01_rbdry);
[deri_plus_Pi11, deri_minus_Pi11] = deri(Pi11, Pi11_lbdry, Pi11_rbdry);
[deri_plus_Phi00, deri_minus_Phi00] = deri(Phi00, Phi00_lbdry, Phi00_rbdry);
[deri_plus_Phi01, deri_minus_Phi01] = deri(Phi01, Phi01_lbdry, Phi01_rbdry);
[deri_plus_Phi11, deri_minus_Phi11] = deri(Phi11, Phi11_lbdry, Phi11_rbdry);
if (USE_S)
    [deri_plus_S, deri_minus_S] = deri(S, S_lbdry, S_rbdry);
    [deri_plus_Pi_S, deri_minus_Pi_S] = deri(Pi_S, Pi_S_lbdry, Pi_S_rbdry);
    [deri_plus_Phi_S, deri_minus_Phi_S] = deri(Phi_S, Phi_S_lbdry, Phi_S_rbdry);
end;
if (USE_psi)
    [deri_plus_psi, deri_minus_psi] = deri(psi, psi_lbdry, psi_rbdry);
    [deri_plus_Pi_psi, deri_minus_Pi_psi] = deri(Pi_psi, Pi_psi_lbdry, Pi_psi_rbdry);
    [deri_plus_Phi_psi, deri_minus_Phi_psi] = deri(Phi_psi, Phi_psi_lbdry, Phi_psi_rbdry);
end;

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
if (USE_S)
    avg_deri_S = 0.5.*(deri_plus_S+deri_minus_S);
    avg_deri_Pi_S = 0.5.*(deri_plus_Pi_S+deri_minus_Pi_S);
    avg_deri_Phi_S = 0.5.*(deri_plus_Phi_S+deri_minus_Phi_S);
end;
if (USE_psi)
    avg_deri_psi = 0.5.*(deri_plus_psi+deri_minus_psi);
    avg_deri_Pi_psi = 0.5.*(deri_plus_Pi_psi+deri_minus_Pi_psi);
    avg_deri_Phi_psi = 0.5.*(deri_plus_Phi_psi+deri_minus_Phi_psi);
end;

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
if (USE_S)
    Hhat_S = -(1+paragamma1).*shift.*avg_deri_S - max_alpha.*0.5.*(deri_plus_S-deri_minus_S);
    Hhat_Pi_S = -paragamma1.*paragamma2.*shift.*avg_deri_S - shift.*avg_deri_Pi_S + lapse.*gamma11.*avg_deri_Phi_S ...
             - max_alpha.*0.5.*(deri_plus_Pi_S-deri_minus_Pi_S);
    Hhat_Phi_S = -paragamma2.*lapse.*avg_deri_S + lapse.*avg_deri_Pi_S - shift.*avg_deri_Phi_S ...
             - max_alpha.*0.5.*(deri_plus_Phi_S-deri_minus_Phi_S);
end;
if (USE_psi)
    Hhat_psi = -(1+paragamma1).*shift.*avg_deri_psi - max_alpha.*0.5.*(deri_plus_psi-deri_minus_psi);
    Hhat_Pi_psi = -paragamma1.*paragamma2.*shift.*avg_deri_psi - shift.*avg_deri_Pi_psi ...
            + lapse.*gamma11.*avg_deri_Phi_psi - max_alpha.*0.5.*(deri_plus_Pi_psi-deri_minus_Pi_psi);
    Hhat_Phi_psi = -paragamma2.*lapse.*avg_deri_psi + lapse.*avg_deri_Pi_psi - shift.*avg_deri_Phi_psi ...
            - max_alpha.*0.5.*(deri_plus_Phi_psi-deri_minus_Phi_psi);
end;

%RHS
rhs_g00 = src_g00 - Hhat_g00;
rhs_g01 = src_g01 - Hhat_g01;
rhs_g11 = src_g11 - Hhat_g11;
rhs_Pi00 = src_Pi00 - Hhat_Pi00;
rhs_Pi01 = src_Pi01 - Hhat_Pi01;
rhs_Pi11 = src_Pi11 - Hhat_Pi11;
rhs_Phi00 = src_Phi00 - Hhat_Phi00;
rhs_Phi01 = src_Phi01 - Hhat_Phi01;
rhs_Phi11 = 0.*x; %src_Phi11 - Hhat_Phi11;
if (USE_S)
    rhs_S = src_S - Hhat_S;
    rhs_Pi_S = src_Pi_S - Hhat_Pi_S;
    rhs_Phi_S = src_Phi_S - Hhat_Phi_S;
end;
if (USE_psi)
    rhs_psi = src_psi - Hhat_psi;
    rhs_Pi_psi = src_Pi_psi - Hhat_Pi_psi;
    rhs_Phi_psi = src_Phi_psi - Hhat_Phi_psi;
end;

%constraints
Cr00 = avg_deri_g00 - Phi00;
Cr01 = avg_deri_g01 - Phi01;
Cr11 = avg_deri_g11 - Phi11;
return
