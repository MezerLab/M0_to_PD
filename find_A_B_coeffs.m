function [A,B] = find_A_B_coeffs(PD,map,BM)


mdl=fitlm(1./PD(BM),map(BM));
A(1) = mdl.Coefficients.Estimate(2);
B(1) = mdl.Coefficients.Estimate(1);