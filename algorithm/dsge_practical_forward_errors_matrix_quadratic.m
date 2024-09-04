function [results] = dsge_practical_forward_errors_matrix_quadratic(inputs)
%Matrix quadratic backward errors and conditioning number follow Higham and Kim (2001)
%0=A*X^2+B*X+C
%X is the n by n solvent
%A, B, and C are n by n matrices
%Alexander Meyer-Gohde
%24/03/2022
if inputs.nstatic>0
A=[inputs.A_static;inputs.AA];%sparse(inputs.A_full);
B=[inputs.B_static inputs.B_rest;zeros(inputs.ndynamic,inputs.nstatic) inputs.BB];%sparse(inputs.B_full);
C=[inputs.C_static;inputs.CC];%sparse(inputs.C_full);
else
    A=inputs.AA;
    B=inputs.BB;
    C=inputs.CC;
end
A=[zeros(inputs.endo_nbr,inputs.npred+inputs.nstatic) A];
C=[zeros(inputs.endo_nbr,inputs.nstatic) C zeros(inputs.endo_nbr,inputs.nfwrd)];

P=inputs.X;
F=A*P+B;
results=NaN(8,5);
P_F=norm(P,'fro');
R_P=F*P+C;
R_P_F=norm(R_P,'fro');
[ny,~]=size(P);
[temp_P_FE in difp]=SYLG_allocated(ny,ny,real(F),eye(ny),A,real(P)',real(R_P),1);
if ny<51
V=kron(eye(ny),F)+kron(P',A);
results(4,1)=svds(V,1,'smallest');
else
    results(4,1)=difp;
end
results(7,1)=norm(temp_P_FE)/P_F;
results(8,1)=1/results(4,1)*(R_P_F/P_F);




