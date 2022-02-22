function apply_freezing
Globals1D;
GlobalsGR;

%inner boundary
%rhs_g00_in = rhs_g00(vmapI);
%rhs_g01_in = rhs_g01(vmapI);
%rhs_g11_in = rhs_g11(vmapI);
% 
%rhs_Pi00_in = rhs_Pi00(vmapI);
%rhs_Pi01_in = rhs_Pi01(vmapI);
%rhs_Pi11_in = rhs_Pi11(vmapI);
%                
%rhs_Phi00_in = rhs_Phi00(vmapI);
%rhs_Phi01_in = rhs_Phi01(vmapI);
%rhs_Phi11_in = rhs_Phi11(vmapI);
%
%rhs_g00(vmapI) = 0.;
%rhs_g01(vmapI) = 0.;
%rhs_g11(vmapI) = 0.;
%
%rhs_Pi00(vmapI) = rhs_Pi00_in - paragamma2*rhs_g00_in;
%rhs_Pi01(vmapI) = rhs_Pi01_in - paragamma2*rhs_g01_in;
%rhs_Pi11(vmapI) = rhs_Pi11_in - paragamma2*rhs_g11_in;
  
%rhs_Phi00(vmapI);
%rhs_Phi01(vmapI);
%rhs_Phi11(vmapI);

%out bdry
rhs_g00_out = rhs_g00(vmapO);
rhs_g01_out = rhs_g01(vmapO);
rhs_g11_out = rhs_g11(vmapO);
 
rhs_Pi00_out = rhs_Pi00(vmapO);
rhs_Pi01_out = rhs_Pi01(vmapO);
rhs_Pi11_out = rhs_Pi11(vmapO);
                
rhs_Phi00_out = rhs_Phi00(vmapO);
rhs_Phi01_out = rhs_Phi01(vmapO);
rhs_Phi11_out = rhs_Phi11(vmapO);

%rhs_S_out = rhs_S(vmapO);
%rhs_Pi_S_out = rhs_Pi_S(vmapO);
%rhs_Phi_S_out = rhs_Phi_S(vmapO);

b = sqrt(1/g11(vmapO));
%rhs_g00(vmapO) = 0.;
%rhs_g01(vmapO) = 0.;
%rhs_g11(vmapO) = 0.;

rhs_Pi00(vmapO) = paragamma2/2*rhs_g00_out + 1/2*rhs_Pi00_out ...
        + b/2*rhs_Phi00_out;% - 0.25*(Phi00(vmapO)-Phi00_exact(vmapO))*power(b,3)*rhs_g11_out;
rhs_Pi01(vmapO) = paragamma2/2*rhs_g01_out + 1/2*rhs_Pi01_out ...
        + b/2*rhs_Phi01_out;% - 0.25*(Phi01(vmapO)-Phi01_exact(vmapO))*power(b,3)*rhs_g11_out;
rhs_Pi11(vmapO) = paragamma2/2*rhs_g11_out + 1/2*rhs_Pi11_out ...
        + b/2*rhs_Phi11_out;% - 0.25*(Phi11(vmapO)-Phi11_exact(vmapO))*power(b,3)*rhs_g11_out;
rhs_Phi00(vmapO) = -paragamma2/2/b*rhs_g00_out ...
        + 1/2/b*rhs_Pi00_out + 1/2*rhs_Phi00_out;% + 0.25*(Phi00(vmapO)-Phi00_exact(vmapO))*power(b,2)*rhs_g11_out;
rhs_Phi01(vmapO) = -paragamma2/2/b*rhs_g01_out ...
        + 1/2/b*rhs_Pi01_out + 1/2*rhs_Phi01_out;% + 0.25*(Phi01(vmapO)-Phi01_exact(vmapO))*power(b,2)*rhs_g11_out;
rhs_Phi11(vmapO) = -paragamma2/2/b*rhs_g11_out ...
        + 1/2/b*rhs_Pi11_out + 1/2*rhs_Phi11_out;% + 0.25*(Phi11(vmapO)-Phi11_exact(vmapO))*power(b,2)*rhs_g11_out;

%rhs_S(vmapO) = 0.;
%rhs_Pi_S(vmapO) = -paragamma2/2*rhs_S_out + 1/2*rhs_Pi_S_out ...
%        + b/2*rhs_Phi_S_out ;
%rhs_Pi_psi(vmapO) = paragamma2/2*rhs_psi_out + 1/2*rhs_Pi_psi_out ...
%        + b/2*rhs_Phi_psi_out ;
%
%rhs_Phi_S(vmapO) = -paragamma2/2/b*rhs_S_out ...
%        + 1/2/b*rhs_Pi_S_out + 1/2*rhs_Phi_S_out ;
%rhs_Phi_psi(vmapO) = -paragamma2/2/b*rhs_psi_out ...
%        + 1/2/b*rhs_Pi_psi_out + 1/2*rhs_Phi_psi_out ;

%rhs_Pi00(vmapO) = 0.;
%rhs_Pi01(vmapO) = 0.; 
%rhs_Pi11(vmapO) = 0.; 
%rhs_Phi00(vmapO) = 0.;
%rhs_Phi01(vmapO) = 0.;
%rhs_Phi11(vmapO) = 0.;

%rhs_S = 0.*x;
%rhs_Pi_S = 0.*x;
%rhs_Phi_S = 0.*x;
%rhs_psi = 0.*x;
%rhs_Pi_psi = 0.*x;
%rhs_Phi_psi = 0.*x;
