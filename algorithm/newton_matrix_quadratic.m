function [X,varargout] = newton_matrix_quadratic(matrix_quadratic,varargin)
% A returns a solvent of the the matrix quadratic equation
% A X^2 +B X +C=0
%
% INPUTS
%   matrix_quadratic    [structure] A structure containing the square
%                                   matrices A, B, and C
%
%   options (optional)  [structure] A structure containing the following
%                                   possible options:
%
% OUTPUTS
%   X                   [matrix]    A solvent of 0=A*X^2+B*X+C
%
%   output (optional)   [structure] A structure containing the following
%                                   possible outputs:
%
%   diff (optional)     [scalar]    The last value of the convergence
%                                   criteron
%
%   j    (optional)     [scalar]    The number of iterations
%
%   resid (optional)    [matrix]    The residual of A*X^2+B*X+C
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

% Initialization. Housekeeping: Extract the matrices from the struct
A=matrix_quadratic.AA;
B=matrix_quadratic.BB;
C=matrix_quadratic.CC;

nfwrd=matrix_quadratic.nfwrd;
npred=matrix_quadratic.npred;
nboth=matrix_quadratic.nboth;
nsfwrd=matrix_quadratic.nsfwrd;
nspred=matrix_quadratic.nspred;
ndynamic=matrix_quadratic.ndynamic;

% % Check that the matrices A, B, and C are square and of the same dimension
% [n,a2]=size(A);
% [b1,b2]=size(B);
% [c1,c2]=size(C);
% 
% if n~=a2
%     disp('A is not a square matrix')
%     return;
% end
% if b1~=b2
%     disp('B is not a square matrix')
%     return;
% end
% if c1~=c2
%     disp('C is not a square matrix')
%     return;
% end
% if n~=b1|| b1~=c1
%     disp('A, B, and C are not of the same dimension')
%     return;
% end

% Check the number of input arguments
if nargin>2
    disp('Too many input arguments')
    return;
else    % allocate the correct options depending on the chosen algorithm
    if nargin==2
        options=varargin{1};
    end
    % 'baseline' is the default algorithm
    if nargin==2 && isfield(options,"algorithm")
        algorithm=options.algorithm;
    else
        algorithm="baseline";
    end
    % determine the sylvester method 
    if nargin==2 && isfield(options,"sylvester_method")
        sylvester_method=options.sylvester_method;
    else
        sylvester_method="hessenberg-schur";
    end
    % determine the convergence metric 
    if nargin==2 && isfield(options,"convergence_metric")
        convergence_metric=options.convergence_metric;
    else
        convergence_metric="residual";
    end
    % determine if dynare's sylvester solver is used
    if nargin==2 && isfield(options,"dynare_reduced_sylvester")
        dynare_reduced_sylvester=options.dynare_reduced_sylvester;
    else
        dynare_reduced_sylvester=0;
    end
    % determine the initial guess of the Pj solution matrix
    if nargin==2 && isfield(options,"initial")
        Pj=options.initial;
    else
        Pj=zeros(ndynamic,nspred);
    end
    % determine the maximum amount of iterations
    if nargin==2 && isfield(options,"maximum_iterations")
        maxiter=options.maximum_iterations;
    else
        maxiter=100;
    end
    % determine the convergence tolerance
    if nargin==2 && isfield(options,"convergence_tolerance")
        tol=options.convergence_tolerance;
    else
        tol=ndynamic*eps;%Default convergence criterion, p. 508 Bai, Z.-Z. and Gao,
        %Y.-H. (2007). Modified bernoulli iteration methods for
        %quadratic matrix equation. Journal of Computational
        %Mathematics, 25(5):498�511
    end
end

%no coefficient update in algorithm 'modified', call initial matrix guess
%Pmod
if strcmp(algorithm,'modified')
    coefficient_update=0;
    if nargin==2 && isfield(options,"initial")
        Pmod=options.initial;
    else
        Pmod=zeros(ndynamic,nspred);
    end
else
    coefficient_update=1;
end

%determine if the algorithm uses line searches and choose line search
%method
if strcmp(algorithm,'line_search')||strcmp(algorithm,'occ_line_search')||strcmp(algorithm,'occ_line_search_samanskii')
    line_search=1;
    if nargin==2 && isfield(options,"line_search_method")
        line_search_method=options.line_search_method;
    else
        line_search_method='direct';
    end
else
    line_search=0;
end

%determine if the algorithm uses occasional line searches and choose
%convergence metric
if strcmp(algorithm,'occ_line_search')||strcmp(algorithm,'occ_line_search_samanskii')
    occasional_line_search=1;
    if nargin==2 && isfield(options,"convergence_metric_occ")
        convergence_metric_occ=options.convergence_metric_occ;
    else
        convergence_metric_occ='residual';
    end
    if nargin==2 && isfield(options,"tol_occ")
        tol_occ=options.tol_occ;
    else
        tol_occ=tol;
    end
else
    occasional_line_search=0;
end

%determine if a samanskii-step is made
if strcmp(algorithm,'samanskii')||strcmp(algorithm,'occ_line_search_samanskii')
    samanskii=1;
    if nargin==2 && isfield(options,"samanskii_steps")
        samanskii_steps=options.samanskii_steps;
    else
        samanskii_steps=1;
    end
else
    samanskii=0;
end


% %reorder matrices if using Dynare reduced sylvester solvers
% if dynare_reduced_sylvester
%     cols_with_all_zeros = find(all(A==0));
%     cols_with_not_all_zeros = find(any(A~=0));
%     change_order=[cols_with_all_zeros cols_with_not_all_zeros];
%     Ar=A(change_order,cols_with_not_all_zeros);
%     A=A(change_order,change_order);
%     B=B(change_order,change_order);
%     C=C(change_order,change_order);
%     Pj=Pj(change_order,change_order);
%     reverse_order=1:n;reverse_order(change_order)=reverse_order;
%     if ~coefficient_update
%         Pmod=Pmod(change_order,change_order);
%     end
% end


% Initialization
diff=tol+1;                     % Initializing convergence check
j=0;                            % Counter for number of iterations
G=[A*Pj(npred+1:end,:) zeros(ndynamic,nfwrd)]+B;                  		% call the bracket G
M=G*Pj+C;                       % calculate M 
P0=Pj;                          
tj=1;                           % default Newton step size is 1

%Begin algorithm(s)
while diff>tol
    % determine if matrix Pmod is updated or not based on algorithm
    if coefficient_update
        Pmod=Pj;
    end
    % solve for matrix dPj depending on chosen solution method
    [U,T]=qr(G(:,1:npred));
T=T(1:size(T,2),:);
UT=U';
UT_1= UT(1:size(T,2),:);
UT_2= UT(size(T,2)+1:end,:);
UT2_G=UT_2*A;
UT2_B=[UT2_G*Pmod(npred+1:end,npred+1:npred+nboth) zeros(nsfwrd,nfwrd)]+UT_2*B(:,npred+1:end);
UT2_T=UT_2*M(:,1:npred+nboth);
[dPj]=generalized_sylvester_solvers(UT2_B, UT2_G, Pmod(1:npred+nboth,1:npred+nboth), UT2_T,sylvester_method);

T1_U1=T\UT_1;
TU1_G=T1_U1*A;
TU1_B=[TU1_G*Pmod(npred+1:end,npred+1:npred+nboth) zeros(npred,nfwrd)]+T1_U1*B(:,npred+1:end);
TU1_T=T1_U1*M(:,1:npred+nboth);
dPj=[-TU1_T-TU1_B*dPj-TU1_G*dPj*Pmod(1:npred+nboth,1:npred+nboth);dPj];
%dPj=[real(dPj) zeros(ndynamic,nfwrd)];
 %   [dPj]=generalized_sylvester_solvers(G, A, Pmod, M,sylvester_method);
    if occasional_line_search
        % determine if line searches are done based on a convergence
        % criterion 
        Pjt=Pj+dPj;
        Gt=[A*Pjt(npred+1:end,:) zeros(ndynamic,nfwrd)]+B;
        Mt=Gt*Pjt+C;
        if convergence_criterion(Mt,Pjt,Pj,Gt,A,B,C,convergence_metric_occ)>tol_occ
            line_search=1;
        else
            line_search=0;
        end
    end
    if line_search                  % use chosen line search method to determine scalar tj 
                                    % for applicable algorithms
        alpha_hk=norm(M,'fro')^2;
        gamma_hk=dPj(npred+1:end,:)*dPj(1:nspred,:);
        A_dP_j_2=A*gamma_hk;
        if matrix_quadratic.nstatic>0
            gamma_hk=matrix_quadratic.A_static*gamma_hk;
            gamma_hk=trace(gamma_hk*gamma_hk')+trace(A_dP_j_2*A_dP_j_2');
        else
            gamma_hk=trace(A_dP_j_2*A_dP_j_2');
        end
        beta_hk=M'*A_dP_j_2;
        beta_hk=2*trace(beta_hk);
        [tj] = line_search_solvers(alpha_hk,beta_hk,gamma_hk,line_search_method);
% Pj_long=[Pmod zeros(ndynamic,nfwrd)];
% matrix_quadratic.X=Pj_long; [matrix_quadratic]=complete_reduced_matrix_quadratic(matrix_quadratic);
% M_full=(matrix_quadratic.A_full*matrix_quadratic.X+matrix_quadratic.B_full)*matrix_quadratic.X+matrix_quadratic.C_full;
% alpha_hk_full=norm(M_full,'fro')^2;
% dPj_long=[dPj zeros(ndynamic,nfwrd)];
% matrix_quadratic.X=dPj_long; [matrix_quadratic]=complete_reduced_matrix_quadratic(matrix_quadratic);
% A_dP_j_2_full=matrix_quadratic.A_full*matrix_quadratic.X*matrix_quadratic.X;
% gamma_hk_full=norm(A_dP_j_2_full,'fro')^2;
% beta_hk_full=M_full'*A_dP_j_2_full;
% beta_hk_full=2*trace(beta_hk_full);
% [tj] = line_search_solvers_deflate(alpha_hk_full,beta_hk_full,gamma_hk_full,line_search_method);

    end
    Pj=P0+tj*dPj;                   % calculate solution matrix Pj
    %M=(A*Pj+B)*Pj+C;                % use solution to update matrix M
    M=([A*Pj(npred+1:end,:) zeros(ndynamic,nfwrd)]+B)*Pj+C;  
    if samanskii && ~line_search
        samanskii_step=1;
%        tj=1;
%        Pmod=P0;
        while samanskii_step<=samanskii_steps

UT2_T=UT_2*M(:,1:npred+nboth);
[dPj]=generalized_sylvester_solvers(UT2_B, UT2_G, Pmod(1:npred+nboth,1:npred+nboth), UT2_T,sylvester_method);
TU1_T=T1_U1*M(:,1:npred+nboth);
dPj=[-TU1_T-TU1_B*dPj-TU1_G*dPj*Pmod(1:npred+nboth,1:npred+nboth);dPj];
%dPj=[real(dPj) zeros(ndynamic,nfwrd)];
            Pj=Pj+dPj;                   % calculate solution matrix Pj
            M=([A*Pj(npred+1:end,:) zeros(ndynamic,nfwrd)]+B)*Pj+C;             % use solution to update matrix M
            samanskii_step=samanskii_step+1;
        end
    end
    % update coefficient G for relevant algorithms
    if coefficient_update          
        G=[A*Pj(npred+1:end,:) zeros(ndynamic,nfwrd)]+B;                        
    end
    %     if close_to_convergence_behavior
    %         %criterion
    %         %H&K turn off line searches when close to convergence...
    %     end
    diff=convergence_criterion(M,Pj,P0,G,A,B,C,convergence_metric);           % determine convergence criterion
    j=j+1;                               % advance the counter
    P0=Pj;
    if j>maxiter
        disp('Maximum iterations reached')
        break
    end
end
%End algorithm(s)

% reorder solution matrix Pj again when using dynare routines and call
% output X
% if dynare_reduced_sylvester
%     X=Pj(reverse_order,reverse_order);
% else
    X=real([Pj zeros(ndynamic,nfwrd)]);
%end

% determine optional output
if nargout>2
    disp('Too many output arguments')
    return;
elseif nargout ==2
    %output.diff=diff;
    %output.j=j;
    %output.M=([A*Pj(1:nsfwrd,:) zeros(ndynamic,nfwrd)]+B)*Pj+C;
    %varargout{1}=output;
    varargout{1}=j;
end

end

