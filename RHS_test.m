function [rhs_A, rhs_B] = RHS_test(A, B, A_bdry, B_bdry)
Globals1D;

M_cos = tanh(x-2);

src_A = (-1/2+M_cos./2).*(-2).*power(B,3) + (1/2+M_cos./2).*(-1).*power(x,-2);
src_B = (1/2+M_cos./2).*(-2).*power(B,3) + (-1/2+M_cos./2).*(-1).*power(x,-2);

%Hhat terms
[deri_plus_A, deri_minus_A] = deri(A, A, A_bdry);
[deri_plus_B, deri_minus_B] = deri(B, B, B_bdry);

max_alpha = 1; 
avg_deri_A = 0.5.*(deri_plus_A+deri_minus_A);
avg_deri_B = 0.5.*(deri_plus_B+deri_minus_B);

Hhat_A = (-1/2+M_cos./2).*avg_deri_A + (1/2+M_cos./2).*avg_deri_B - max_alpha.*0.5.*(deri_plus_A-deri_minus_A);
Hhat_B = (1/2+M_cos./2).*avg_deri_A + (-1/2+M_cos./2).*avg_deri_B - max_alpha.*0.5.*(deri_plus_B-deri_minus_B);

%RHS
rhs_A = src_A - Hhat_A;
rhs_B = src_B - Hhat_B;

return
