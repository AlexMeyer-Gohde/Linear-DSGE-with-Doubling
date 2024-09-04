function [Y]=generalized_sylvester_solvers(A,B,C,D,method) 
% Y returns a solution to the equation A * Y + B * Y * C + D=0
%
% INPUTS
%   A, B, C, D  [matrix]     matrices A, B, C, D from the above equation                                  
%
%   method      [string]     a string determining which sylvester solver
%                            is used
%                                   
%
% OUTPUTS
%   Y                   [matrix]    A solvent of A * Y + B * Y * C + D=0
%
%
% ALGORITHMS
%   Meyer-Gohde, Alexander and Saecker, Johanna (2022). SOLVING LINEAR DSGE
%   MODELS WITH NEWTON METHODS
%
%
% Copyright (C) 2022 Alexander Meyer-Gohde & Johanna Saecker
%
% This is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.


% the 'hessenberg-schur' option follows Higham & Kim, 2001 (equations below
% refer to their paper)
n1=length(A);
n2=length(C);
if strcmp(method,'hessenberg-schur')
    Y=zeros(n1,n2);                   % create Y-matrix
    TY=zeros(n1,n2);
    
    [H,T,Ws,Z]=hess(A,B);           % equation (3.3) with Hessenberg-triangular H instead of S
    [U,R]=schur(C,'complex');                 % equation (3.3)
    F=-Ws*D*U;                      % equation (3.4)
    for k=1:n2                       % loop through columns of Y
        lhs=R(k,k)*T;
        lhs=lhs+H;
        rhs=TY(:,1:k-1)*R(1:k-1,k);
        rhs=F(:,k)-rhs;
        Y(:,k)=lhs\rhs;
        TY(:,k)=T*Y(:,k);
    end
    %  max(max(abs(Y_amg-Y)))
    Y=Z*Y*U';

% the 'dynare' option uses dynare's generalized sylvester solver, solving
% the equation A_ * Y + B_ * Y * C_ + D_ = 0
% caution: for dynare versions older than 5.0 outputs seem to be inverted
elseif strcmp(method,'dynare')
    [Y,err] = gensylv(1, A, B, C, -D); %Dynare [err, ghx_other] = gensylv(1, A_, B_, C_,-D_); solves  A_ * x + B_ * x * C_ + D_ = 0

% the 'dlyap' option uses Matlab's solver for discrete-time Lyapunov
% equations solving equation A_t * Y * C - Y + D_t = 0
elseif strcmp(method,'dlyap')
    A_t=-A\B; 
    D_t=-A\D; 
    Y=dlyap(A_t,C,D_t);

elseif strcmp(method,'dlyap_stripped')    
    
    CC=-A\D;
    A=-A\B;
    B=C;
   % X = dsylvester(A, B, -C)
 %  X = dsylvester(A, B, C)
% Solve Sylvester Equation A*X*B - I = C.
  % [sA,pA,A] = mscale(A,'fullbal');
  % [sB,pB,B] = mscale(B,'fullbal');
  % CC(pA,pB) = lrscale(CC,1./sA,sB);  % TA\C*TB


        [QA, TA] = schur(A);
        CC = QA'*CC;
        [QB, TB] = schur(B);
        CC = CC*QB;

    
    % Solve Sylvester Equation -TA*X*TB' + X = -QA'*C*QB.
    Y = matlab.internal.math.sylvester_tri(TA, 'I', CC, 'I', -TB, 'notransp');
 Y = QA*Y;
 Y = Y*QB';
 %Y = lrscale(Y(pA,pB),sA,1./sB);  % TA*X/TB

elseif strcmp(method,'lapack')
    % ZTGSYL( TRANS, IJOB, M, N, A, LDA, B, LDB, C, LDC, D,
%      $                   LDD, E, LDE, F, LDF, SCALE, DIF, WORK, LWORK,
%      $                   IWORK, INFO )

%     A * R - L * B = scale * C )
%                               )                                 (1)
%     D * R - L * E = scale * F )
%  where A and D are M-by-M matrices, B and E are N-by-N matrices and
%  C, F, R and L are M-by-N matrices

[A1,D1,P,Q] = qz(A,B);
%[B1,E1,U,V_T] = qz(B,E);
[U,B1] = schur(-C,'complex');

C1=-P*D*U;
F1=zeros(size(A1,1),size(B1,1));
E1=eye(size(B1));
%results = lapack('ZTGSYL', 'N', 0, size(A1,2), size(B1,2), A1, size(A1,1), B1, size(B1,1), C1, size(C1,1), D1, size(D1,1), E1, size(E1,1), F1, size(F1,1), 1, 1, 1, 1, ones(size(A1,2)+size(B1,2)+6,1), 0);
%Y=Q*results{9}*U';
%[R L dif scale work iwork] = ztgsyl('N',0,A1,B1,C1,D1,E1,F1);
[R, ~, ~, ~, ~, ~] = ztgsyl('N',0,complex(A1),complex(B1),complex(C1),complex(D1),complex(E1),complex(F1));
Y=Q*R*U';
elseif strcmp(method,'slicot')

    CC=-A\D;
    A=A\B;
    B=C;
    [N, M]=size(CC);
    [Y,~,~] = SB04QD(N,M,A,B,CC);
Y=Y(1:N,1:M);
else
    disp('Unknown Sylvester Solver');
end
