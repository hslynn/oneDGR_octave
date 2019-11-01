function rhs_A = RHS_testA(A, A_exact)
Globals1D;

M_tanh = tanh(x-2);

src_A = M_tanh.*(-1).*power(x,-2);

%Hhat terms
[deri_plus_A, deri_minus_A] = deri(A, A, A);

max_alpha = max(max(abs(M_tanh)));
avg_deri_A = 0.5.*(deri_plus_A+deri_minus_A);

Hhat_A = M_tanh.*avg_deri_A - max_alpha.*0.5.*(deri_plus_A-deri_minus_A);

%RHS
rhs_A = src_A - Hhat_A;
%rhs_A(0.5*mapO-9) = 0;
%rhs_A(0.5*mapO+10) = 0;
return
