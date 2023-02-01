function [A,B] = find_A_B_coeffs(PD,map,BM)

warning off
mdl=fitlm(1./PD(BM),map(BM));
warning on
A(1) = mdl.Coefficients.Estimate(2);
B(1) = mdl.Coefficients.Estimate(1);
