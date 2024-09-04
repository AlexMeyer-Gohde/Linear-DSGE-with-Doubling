function [X,j] = bernoulli_matrix_quadratic_fast(matrix_quadratic,X,options)
% Returns a (the minimal) solvent X of the matrix quadratic equation
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
%   Meyer-Gohde, Alexander (2022). SOLVING LINEAR DSGE
%   MODELS WITH BERNOULLI METHODS
%
%
% Copyright (C) 2022 Alexander Meyer-Gohde
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


% % Check the number of input arguments
% if nargin>2
%     disp('Too many input arguments')
%     return;
% else    % allocate the correct options depending on the chosen algorithm
%     if nargin==2
%         options=varargin{1};
%     end
%     if nargin==2 && isfield(options,"baseline")
%         baseline=options.baseline;
%     else
%         baseline=1;
%     end
%     if nargin==2 && ~baseline
%         if isfield(options,"mbi")
%             mbi=options.mbi;
%         else
%             mbi=0;
%         end
%         if isfield(options,"newton")
%             newton=1;
%             if strcmp(options.newton,'row')
%                 newton_row=1; newton_column=0; newton_block=0;newton_opt=0;
%             elseif strcmp(options.newton,'column')
%                 newton_column=1; newton_row=0; newton_block=0;newton_opt=0;
%             elseif strcmp(options.newton,'block')
%                 newton_block=1;newton_column=0; newton_row=0;newton_opt=0;
%             else
%                 newton_block=0;newton_column=0; newton_row=0; newton_opt=1;
%             end
%             if newton && isfield(options,"sylvester_method")
%                 sylvester_method=options.sylvester_method;
%             else
%                 sylvester_method="dlyap_stripped";
%             end
%         else
%             newton=0;
%         end
%     end
%     if nargin==2 && isfield(options,"line_search")
%         line_search=options.line_search;
%     else
%         line_search=0;
%     end
%     if nargin==2 && isfield(options,"convergence_metric")
%         convergence_metric=options.convergence_metric;
%         if strcmp(convergence_metric,'reldiff')
%             convergence_metric=1;
%         elseif strcmp(convergence_metric,'residual')
%             convergence_metric=2;
%         elseif strcmp(convergence_metric,'fe1')
%             convergence_metric=3;
%         elseif strcmp(convergence_metric,'fe1_sparse')
%             convergence_metric=4;
%         end
%     else
%         convergence_metric=1;
%     end
%     % determine the initial guess of the Pj solution matrix
%     if nargin==2 && isfield(options,"initial")
%         X=options.initial;
%     else
%         X=zeros(ndynamic,nspred);
%     end
%     % determine the maximum amount of iterations
%     if nargin==2 && isfield(options,"maximum_iterations")
%         maxiter=options.maximum_iterations;
%     else
%         maxiter=100;
%     end
%     % determine the convergence tolerance
%     if nargin==2 && isfield(options,"convergence_tolerance")
%         tol=options.convergence_tolerance;
%     else
%         tol=ndynamic*eps;%Default convergence criterion
%     end
%     if nargin==2 && isfield(options,"maximum_restarts")
%         max_restart=options.maximum_restarts;
%     else
%         max_restart=1;
%     end
%     if nargin==2 && isfield(options,"newton_power")
%         newton_power=options.newton_power;
%     else
%         newton_power=1;
%     end
% end

% baseline=0;
% mbi=0;
% newton=0;
% newton_block=0;
% newton_comlumn=0;
% newton_row=0;
% newton_opt=0;
% newton_power=1;
% line_search=1;
% convergence_metric=4;
% maxiter=100;
% tol=ndynamic*eps;
% max_restart=0;

baseline=options(1);
mbi=options(2);
newton=options(3);
newton_block=options(4);
newton_column=options(5);
newton_row=options(6);
newton_opt=options(7);
newton_power=options(8);
line_search=options(9);
convergence_metric=options(10);
maxiter=options(11);
tol=options(12);
max_restart=options(13);

sylvester_method="dlyap_stripped";


restart=0;
diff=tol+1; % Initializing convergence check
nfwrd_select=nspred+1:ndynamic;
nsfwrd_select=npred+1:ndynamic;
nspred_select=1:nspred;
B_1=B(:,nspred_select);
B_2=B(:,nfwrd_select);
X_nsfwrd=X(nsfwrd_select,:);AXB_1=B_1+A*X_nsfwrd; AX_B=[AXB_1 B_2];
%AX_B=[A*X(npred+1:end,:) zeros(ndynamic,nfwrd)]+B;
M=AX_B*X+C;
j=0; % Counter for number of iterations
X_0=X;

try
    while diff>tol % As long as the iterations haven't converged...
        if baseline
            X=-AX_B\C;
        elseif mbi
            dlhs = decomposition(AX_B);
            if isIllConditioned(dlhs)
                disp('MBI coefficient breakdown')
                j=1/0;
                break
            else
                X=-dlhs\C; %Functional iteration (Eq. 26, p. 512 from the header)
                for i=2:nspred
                    u=A*(X(npred+1:end,i-1)-X_0(npred+1:end,i-1));
                    p=dlhs\u;
                    eq=X(i-1,i);
                    ep=p(i-1);
                    factor=eq/(1+ep)*p;
                    X(:,i:end)=X(:,i:end)+repmat(factor,[1,nspred+1-i]);
                end
            end
        else
            delta_X=-(AX_B)\(M);
            t=1;
            if line_search
                alpha=sum(sum(M.*M));
                delta_X_nsfwrd=delta_X(nsfwrd_select,:);
                delta_X_nspred=delta_X(nspred_select,:);
                X_nspred=X(nspred_select,:);
                Y_1=delta_X_nsfwrd*delta_X_nspred;
                Z_1=delta_X_nsfwrd*X_nspred;
                Y_2=A*Y_1;
                Z_2=A*Z_1;
                if matrix_quadratic.nstatic>0
                    Y_1=matrix_quadratic.A_static*Y_1;
                    Z_1=matrix_quadratic.A_static*Z_1;
                    gamma=sum(sum(Y_1.*Y_1))+sum(sum(Y_2.*Y_2));
                    beta=sum(sum(Z_1.*Z_1))+sum(sum(Z_2.*Z_2));
                    xi=2*(sum(sum(Y_1.*Z_1))+sum(sum(Y_2.*Z_2)));
                else
                    gamma=sum(sum(Y_2.*Y_2));
                    beta=sum(sum(Z_2.*Z_2));
                    xi=2*sum(sum(Y_2.*Z_2));
                end
                sigma=2*sum(sum(M.*Y_2));
                delta=2*sum(sum(M.*Z_2));
                a=4*gamma;
                b=3*(xi-sigma);
                c=2*(alpha+beta-delta+sigma);
                d=delta-2*alpha;
                p=(3*a*c-b^2)/(3*a^2);
                q=(2*b^3-9*a*b*c+27*a^2*d)/(27*a^3);
                discr=-4*p^3-27*q^2;
                if discr>0+eps %three real roots
                    if (a+b+c+d)*a>0 && a*(3*a+b)>0%No root or two roots in interval 1, inf % No root in interval 1, inf
                        t=1;
                    else
                        k=0:2;
                        root=2*(-p/3)^(1/2)*cos(1/3*acos(3/2*q/p*(-3/p)^(1/2))-2*pi*k/3);
                        root=root-b/(3*a);
                        root=root(root>1);
                        if isempty(root)
                            root=1;
                        end
                        values=(gamma*root.^4+(xi-sigma)*root.^3+(alpha+beta-delta+sigma)*root.^2+(delta-2*alpha)*root+alpha);
                        t=root(values==min(values));
                    end
                elseif discr<0-eps %one real root
                    if (a+b+c+d)*a>0 %No root in interval 1, inf
                        t=1;
                    else
                        if p>0
                            t=-2*(p/3)^(1/2)*sinh(1/3*asinh(3/2*q/p*(3/p)^(1/2)));
                        else
                            t=-2*abs(q)/q*(-p/3)^(1/2)*cosh(1/3*acosh(-3*abs(q)/(2*p)*(-3/p)^(1/2)));
                        end
                        t=t-b/(3*a);
                    end
                else
                    t=-b/(3*a);
                end
                value_t=(gamma*t^4+(xi-sigma)*t^3+(alpha+beta-delta+sigma)*t^2+(delta-2*alpha)*t+alpha);
                value_1=(gamma+xi+beta);
                if (value_1-value_t)/value_t<-eps&&(value_1-value_t)<-eps
                    t=1;
                end
            end
            delta_X=t*delta_X;
            if newton
                tj=1;
                [U,T]=qr(AX_B(:,1:npred));
                T=T(1:size(T,2),:);
                UT=U';
                UT_1= UT(1:size(T,2),:);
                UT_2= UT(size(T,2)+1:end,:);
                UT2_G=UT_2*A;
                UT2_B=[UT2_G*X(npred+1:end,npred+1:npred+nboth) zeros(nsfwrd,nfwrd)]+UT_2*B(:,npred+1:end);
                UT2_T=UT_2*M(:,1:npred+nboth);
                [delta_X_nwt]=generalized_sylvester_solvers(UT2_B, UT2_G, X(1:npred+nboth,1:npred+nboth), UT2_T,sylvester_method);
                T1_U1=T\UT_1;
                TU1_G=T1_U1*A;
                TU1_B=[TU1_G*X(npred+1:end,npred+1:npred+nboth) zeros(npred,nfwrd)]+T1_U1*B(:,npred+1:end);
                TU1_T=T1_U1*M(:,1:npred+nboth);
                delta_X_nwt=[-TU1_T-TU1_B*delta_X_nwt-TU1_G*delta_X_nwt*X(1:npred+nboth,1:npred+nboth);delta_X_nwt];
                if line_search
                    gamma_hk=delta_X_nwt(npred+1:end,:)*delta_X_nwt(1:nspred,:);
                    A_dP_j_2=A*gamma_hk;
                    if matrix_quadratic.nstatic>0
                        gamma_hk=matrix_quadratic.A_static*gamma_hk;
                        gamma_hk=sum(sum(gamma_hk.*gamma_hk))+sum(sum(A_dP_j_2.*A_dP_j_2));
                    else
                        gamma_hk=sum(sum(A_dP_j_2.*A_dP_j_2));
                    end
                    beta_hk=2*sum(sum(M.*A_dP_j_2));
                    [tj] = line_search_solvers(alpha,beta_hk,gamma_hk,'direct');
                end
                delta_X_nwt=tj*delta_X_nwt;
                if newton_block
                    theta=sum(delta_X_nwt(:).*delta_X(:))/(norm(delta_X(:))*norm(delta_X_nwt(:)));
                    theta=max([-1 theta]); theta=min([1 theta]);
                    s=acos(theta)/pi;
                    s=s^newton_power;
                    delta_X=s.*delta_X+(1-s).*delta_X_nwt;
                elseif newton_column
                    EF_sum=sum(delta_X_nwt.*delta_X);
                    theta=EF_sum./(sum(delta_X_nwt.^2).^(1/2).*sum(delta_X.^2).^(1/2));
                    theta(EF_sum==0)=0;
                    theta(theta>1)=1;theta(theta<-1)=-1;
                    s=repmat(acos(theta)/pi,[size(X,1) 1]);
                    s=s.^newton_power;
                    delta_X=s.*delta_X+(1-s).*delta_X_nwt;
                elseif newton_row
                    EF_sum=sum(delta_X_nwt.*delta_X,2);
                    theta=EF_sum./(sum(delta_X_nwt.^2,2).^(1/2).*sum(delta_X.^2,2).^(1/2));
                    theta(EF_sum==0)=0;
                    theta(theta>1)=1;theta(theta<-1)=-1;
                    s=repmat(acos(theta)/pi,[1 size(X,2)]);
                    s=s.^newton_power;
                    delta_X=s.*delta_X+(1-s).*delta_X_nwt;
                elseif newton_opt
                    delta_X_diff=delta_X-delta_X_nwt;
                    a_1=delta_X_nwt(npred+1:end,:)*delta_X_nwt(1:nspred,:);
                    a_2=A*a_1+(1-tj)*M;
                    b_1=delta_X(npred+1:end,:)*X(1:nspred,:)+delta_X_nwt(npred+1:end,:)*delta_X_diff(1:nspred,:)+delta_X_diff(npred+1:end,:)*delta_X_nwt(1:nspred,:);
                    b_2=A*b_1+(tj-t)*M;
                    c_1=delta_X_diff(npred+1:end,:)*delta_X_diff(1:nspred,:);
                    c_2=A*c_1;
                    if matrix_quadratic.nstatic>0
                        a_1=matrix_quadratic.A_static*a_1;
                        b_1=matrix_quadratic.A_static*b_1;
                        c_1=matrix_quadratic.A_static*c_1;
                        a=4*(trace(c_1*c_1')+trace(c_2*c_2'));
                        b=6*(trace(b_1*c_1')+trace(b_2*c_2'));
                        c=2*(trace(a_1*c_1')+trace(a_2*c_2'))+4*(trace(b_1*b_1')+trace(b_2*b_2'));
                        d=2*(trace(a_1*b_1')+trace(a_2*b_2'));
                    else
                        a=4*(trace(c_2*c_2'));
                        b=6*(trace(b_2*c_2'));
                        c=2*(trace(a_2*c_2'))+4*(trace(b_2*b_2'));
                        d=2*(trace(a_2*b_2'));
                    end
                    s=roots([a b c d]);
                    s=s(s>0&s<1);
                    s=[s;0;1];
                    if ~line_search; alpha=norm(M,'fro')^2;end
                    value_s=alpha+s.^4*a/4+s.^3*b/3+s.^2*c/2+s*d;
                    s=s(value_s==min(value_s));
                    s=s(end);

                    delta_X=s.*delta_X+(1-s).*delta_X_nwt;
                end
            end
            X=X+delta_X;
        end
        mustBeNonNan(X)
        X_nsfwrd=X(nsfwrd_select,:);
        AXB_1=B_1+A*X_nsfwrd;
        AX_B=[AXB_1 B_2];
        %AX_B=[A*X(npred+1:end,:) zeros(ndynamic,nfwrd)]+B;
        M=AX_B*X+C;
        %diff=convergence_criterion_fast(M,X,X_0,AX_B,A,B,C,convergence_metric);           % determine convergence criterion
        diff=norm(X-X_0,1)/norm(X,1);
        j=j+1;                               % advance the counter
        X_0=X;
        if j>maxiter
            disp('Maximum iterations reached')
            break
        end
    end


catch
    disp('Encountered numerical difficulties, proceeding more carefully...')
    if rank(AX_B)<ndynamic %If the starting matrix gives a singular coeffcient matrix
        X=lsqminnorm(AX_B,-C);
        AX_B=[A*X(npred+1:end,:) zeros(ndynamic,nfwrd)]+B;
        M=AX_B*X+C;
        j=1; % Counter for number of iterations
        X_0=X;
    end

    while diff>tol % As long as the iterations haven't converged...
        if baseline
            X=-AX_B\C;
        elseif mbi
            dlhs = decomposition(AX_B);
            if isIllConditioned(dlhs)
                disp('MBI coefficient breakdown')
                j=1/0;
                break
            else
                X=-dlhs\C; %Functional iteration (Eq. 26, p. 512 from the header)
                for i=2:nspred
                    u=A*(X(npred+1:end,i-1)-X_0(npred+1:end,i-1));
                    p=dlhs\u;
                    eq=X(i-1,i);
                    ep=p(i-1);
                    factor=eq/(1+ep)*p;
                    X(:,i:end)=X(:,i:end)+repmat(factor,[1,nspred+1-i]);
                end
            end
        else
            delta_X=-(AX_B)\(M);
            t=1;
            if line_search
                alpha=norm(M,'fro')^2;
                delta_X_nsfwrd=delta_X(nsfwrd_select,:);
                delta_X_nspred=delta_X(nspred_select,:);
                X_nspred=X(nspred_select,:);
                Y_1=delta_X_nsfwrd*delta_X_nspred;
                Z_1=delta_X_nsfwrd*X_nspred;
                Y_2=A*Y_1;
                Z_2=A*Z_1;
                if matrix_quadratic.nstatic>0
                    Y_1=matrix_quadratic.A_static*Y_1;
                    Z_1=matrix_quadratic.A_static*Z_1;
                    gamma=sum(sum(Y_1.*Y_1))+sum(sum(Y_2.*Y_2));
                    beta=sum(sum(Z_1.*Z_1))+sum(sum(Z_2.*Z_2));
                    xi=2*(sum(sum(Y_1.*Z_1))+sum(sum(Y_2.*Z_2)));
                else
                    gamma=sum(sum(Y_2.*Y_2));
                    beta=sum(sum(Z_2.*Z_2));
                    xi=2*sum(sum(Y_2.*Z_2));
                end
                sigma=2*sum(sum(M.*Y_2));
                delta=2*sum(sum(M.*Z_2));
                a=4*gamma;
                b=3*(xi-sigma);
                c=2*(alpha+beta-delta+sigma);
                d=delta-2*alpha;
                p=(3*a*c-b^2)/(3*a^2);
                q=(2*b^3-9*a*b*c+27*a^2*d)/(27*a^3);
                discr=-4*p^3-27*q^2;
                if discr>0+eps %three real roots
                    if (a+b+c+d)*a>0 && a*(3*a+b)>0%No root or two roots in interval 1, inf % No root in interval 1, inf
                        t=1;
                    else
                        k=0:2;
                        root=2*(-p/3)^(1/2)*cos(1/3*acos(3/2*q/p*(-3/p)^(1/2))-2*pi*k/3);
                        root=root-b/(3*a);
                        root=root(root>1);
                        if isempty(root)
                            root=1;
                        end
                        values=(gamma*root.^4+(xi-sigma)*root.^3+(alpha+beta-delta+sigma)*root.^2+(delta-2*alpha)*root+alpha);
                        t=root(values==min(values));
                    end
                elseif discr<0-eps %one real root
                    if (a+b+c+d)*a>0 %No root in interval 1, inf
                        t=1;
                    else
                        if p>0
                            t=-2*(p/3)^(1/2)*sinh(1/3*asinh(3/2*q/p*(3/p)^(1/2)));
                        else
                            t=-2*abs(q)/q*(-p/3)^(1/2)*cosh(1/3*acosh(-3*abs(q)/(2*p)*(-3/p)^(1/2)));
                        end
                        t=t-b/(3*a);
                    end
                else
                    t=-b/(3*a);
                end
                value_t=(gamma*t^4+(xi-sigma)*t^3+(alpha+beta-delta+sigma)*t^2+(delta-2*alpha)*t+alpha);
                value_1=(gamma+xi+beta);
                if (value_1-value_t)/value_t<-eps&&(value_1-value_t)<-eps
                    t=1;
                end
            end
            delta_X=t*delta_X;
            if newton
                tj=1;
                [U,T]=qr(AX_B(:,1:npred));
                T=T(1:size(T,2),:);
                UT=U';
                UT_1= UT(1:size(T,2),:);
                UT_2= UT(size(T,2)+1:end,:);
                UT2_G=UT_2*A;
                UT2_B=[UT2_G*X(npred+1:end,npred+1:npred+nboth) zeros(nsfwrd,nfwrd)]+UT_2*B(:,npred+1:end);
                UT2_T=UT_2*M(:,1:npred+nboth);
                [delta_X_nwt]=generalized_sylvester_solvers(UT2_B, UT2_G, X(1:npred+nboth,1:npred+nboth), UT2_T,sylvester_method);
                T1_U1=T\UT_1;
                TU1_G=T1_U1*A;
                TU1_B=[TU1_G*X(npred+1:end,npred+1:npred+nboth) zeros(npred,nfwrd)]+T1_U1*B(:,npred+1:end);
                TU1_T=T1_U1*M(:,1:npred+nboth);
                delta_X_nwt=[-TU1_T-TU1_B*delta_X_nwt-TU1_G*delta_X_nwt*X(1:npred+nboth,1:npred+nboth);delta_X_nwt];
                if line_search
                    gamma_hk=delta_X_nwt(npred+1:end,:)*delta_X_nwt(1:nspred,:);
                    A_dP_j_2=A*gamma_hk;
                    if matrix_quadratic.nstatic>0
                        gamma_hk=matrix_quadratic.A_static*gamma_hk;
                        gamma_hk=sum(sum(gamma_hk.*gamma_hk))+sum(sum(A_dP_j_2.*A_dP_j_2));
                    else
                        gamma_hk=sum(sum(A_dP_j_2.*A_dP_j_2));
                    end
                    beta_hk=2*sum(sum(M.*A_dP_j_2));
                    [tj] = line_search_solvers(alpha,beta_hk,gamma_hk,'direct');
                end
                delta_X_nwt=tj*delta_X_nwt;
                if newton_block
                    theta=sum(delta_X_nwt(:).*delta_X(:))/(norm(delta_X(:))*norm(delta_X_nwt(:)));
                    theta=max([-1 theta]); theta=min([1 theta]);
                    s=acos(theta)/pi;
                    s=s^newton_power;
                    delta_X=s.*delta_X+(1-s).*delta_X_nwt;
                elseif newton_column
                    EF_sum=sum(delta_X_nwt.*delta_X);
                    theta=EF_sum./(sum(delta_X_nwt.^2).^(1/2).*sum(delta_X.^2).^(1/2));
                    theta(EF_sum==0)=0;
                    theta(theta>1)=1;theta(theta<-1)=-1;
                    s=repmat(acos(theta)/pi,[size(X,1) 1]);
                    s=s.^newton_power;
                    delta_X=s.*delta_X+(1-s).*delta_X_nwt;
                elseif newton_row
                    EF_sum=sum(delta_X_nwt.*delta_X,2);
                    theta=EF_sum./(sum(delta_X_nwt.^2,2).^(1/2).*sum(delta_X.^2,2).^(1/2));
                    theta(EF_sum==0)=0;
                    theta(theta>1)=1;theta(theta<-1)=-1;
                    s=repmat(acos(theta)/pi,[1 size(X,2)]);
                    s=s.^newton_power;
                    delta_X=s.*delta_X+(1-s).*delta_X_nwt;
                elseif newton_opt
                    delta_X_diff=delta_X-delta_X_nwt;
                    a_1=delta_X_nwt(npred+1:end,:)*delta_X_nwt(1:nspred,:);
                    a_2=A*a_1+(1-tj)*M;
                    b_1=delta_X(npred+1:end,:)*X(1:nspred,:)+delta_X_nwt(npred+1:end,:)*delta_X_diff(1:nspred,:)+delta_X_diff(npred+1:end,:)*delta_X_nwt(1:nspred,:);
                    b_2=A*b_1+(tj-t)*M;
                    c_1=delta_X_diff(npred+1:end,:)*delta_X_diff(1:nspred,:);
                    c_2=A*c_1;
                    if matrix_quadratic.nstatic>0
                        a_1=matrix_quadratic.A_static*a_1;
                        b_1=matrix_quadratic.A_static*b_1;
                        c_1=matrix_quadratic.A_static*c_1;
                        a=4*(trace(c_1*c_1')+trace(c_2*c_2'));
                        b=6*(trace(b_1*c_1')+trace(b_2*c_2'));
                        c=2*(trace(a_1*c_1')+trace(a_2*c_2'))+4*(trace(b_1*b_1')+trace(b_2*b_2'));
                        d=2*(trace(a_1*b_1')+trace(a_2*b_2'));
                    else
                        a=4*(trace(c_2*c_2'));
                        b=6*(trace(b_2*c_2'));
                        c=2*(trace(a_2*c_2'))+4*(trace(b_2*b_2'));
                        d=2*(trace(a_2*b_2'));
                    end
                    s=roots([a b c d]);
                    s=s(s>0&s<1);
                    s=[s;0;1];
                    if ~line_search; alpha=norm(M,'fro')^2;end
                    value_s=alpha+s.^4*a/4+s.^3*b/3+s.^2*c/2+s*d;
                    s=s(value_s==min(value_s));
                    s=s(end);

                    delta_X=s.*delta_X+(1-s).*delta_X_nwt;
                end
            end
            X=X+delta_X;
        end
        X_nsfwrd=X(nsfwrd_select,:);
        AXB_1=B_1+A*X_nsfwrd;
        AX_B=[AXB_1 B_2];
        %AX_B=[A*X(npred+1:end,:) zeros(ndynamic,nfwrd)]+B;
        M=AX_B*X+C;
        %diff=convergence_criterion_fast(M,X,X_0,AX_B,A,B,C,convergence_metric);           % determine convergence criterion
        diff=norm(X-X_0,1)/norm(X,1);
        j=j+1;                               % advance the counter
        X_0=X;
        if rcond(AX_B)<eps
            restart=restart+1;
            if restart<=max_restart
                X=lsqminnorm(AX_B,-C);
                AX_B=[A*X(npred+1:end,:) zeros(ndynamic,nfwrd)]+B;
                M=AX_B*X+C;
                j=1;
                X_0=X;
            else
                disp('Coefficient breakdown')
                j=1/0;
                break
            end
        end
        if j>maxiter
            disp('Maximum iterations reached')
            break
        end
    end

end


if ~isreal(X)
    X=real(X);
end
X=[X zeros(ndynamic,nfwrd)];
%X(ndynamic,ndynamic)=0;

% % determine optional output
% if nargout>2
%     disp('Too many output arguments')
%     return;
% elseif nargout ==2
%     %output.diff=diff;
%     %output.j=j;
%     %output.M=([A*Pj(1:nsfwrd,:) zeros(ndynamic,nfwrd)]+B)*Pj+C;
%     %varargout{1}=output;
%     varargout=j;
% end




