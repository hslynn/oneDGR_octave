function [rhs_A, rhs_B] = RHS_test_dia(A, B, A_exact, B_exact)
Globals1D;

M_tanh = tanh(x-2);

src_A = (-1).*(-1).*power(B,2);
src_B = M_tanh.*(-1).*power(A,2);

%Hhat terms
[deri_plus_A, deri_minus_A] = deri(A, A(vmapI), A_exact(vmapO));
[deri_plus_B, deri_minus_B] = deri(B, B(vmapI), B_exact(vmapO));

max_alpha = 1.0;
avg_deri_A = 0.5.*(deri_plus_A+deri_minus_A);
avg_deri_B = 0.5.*(deri_plus_B+deri_minus_B);

Hhat_A = (-1).*avg_deri_A - max_alpha.*0.5.*(deri_plus_A-deri_minus_A);
Hhat_B = M_tanh.*avg_deri_B - max_alpha.*0.5.*(deri_plus_B-deri_minus_B);

%RHS
rhs_A = src_A - Hhat_A;
rhs_B = src_B - Hhat_B;

return
