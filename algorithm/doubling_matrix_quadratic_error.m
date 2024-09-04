function [P,varargout] = doubling_matrix_quadratic(matrix_quadratic,varargin)
% Returns a (the minimal) solvent X of the matrix quadratic equation
% A P^2 +B P +C=0
%
% INPUTS
%   matrix_quadratic    [structure] A structure containing the square
%                                   matrices A, B, and C
%
%   options (optional)  [structure] A structure containing options
%
% OUTPUTS
%   P                   [matrix]    A solvent of 0=A*P^2+B*P+C
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
%   Johannes Huber, Alexander Meyer-Gohde, and Johanna Saecker (2023). 
%   SOLVING LINEAR DSGE MODELS WITH STRUCTURE PRESERVING DOUBLING METHODS
%
%
% Copyright (C) 2023 Johannes Huber, Alexander Meyer-Gohde, and 
%   Johanna Saecker
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


% Check the number of input arguments
if nargin>2
    disp('Too many input arguments')
    return;
else    % allocate the correct options depending on the chosen algorithm
    if nargin==2
        options=varargin{1};
    end
end

if nargin==2 && isfield(options,"method")
    if strcmp(options.method,'SF1')
        method=1;
    elseif strcmp(options.method,'SF2')
        method=2;
    else
        disp('Unrecognized method, continuing using SF1');
        method=1;
    end
else
    disp('No method selected, using SF1');
    method=1;
end

if nargin==2 && isfield(options,"P0")
    P0=options.P0;
else
P0=zeros(ndynamic,nspred);
end

if nargin==2 && isfield(options,"convergence_tolerance")
    tol=options.convergence_tolerance;
else
    tol=ndynamic*eps;%Default convergence criterion
end

if nargin==2 && isfield(options,"convergence_metric")
    convergence_metric=options.convergence_metric;
else
    convergence_metric="reldiff";%residual";
end

if nargin==2 && isfield(options,"max_it")
    max_it=options.max_it;
else
    max_it=100;%Default max_it
end
if nargin==2 && isfield(options,"max_restart")
    max_restart=options.max_restart;
else
    max_restart=10;%Default max_restart
end
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
%     else
%         convergence_metric="residual";
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



diff=tol+1;
j=0; 
restart=0;

% initialization
if method==1
    AP=A*P0(npred+1:ndynamic,:);
    AP_B=[AP+B(:,1:nspred),B(:,nspred+1:ndynamic)];
    if rank(AP_B)==ndynamic
        Y=-AP_B\A;
        E=-AP_B\C;
    else
        Y=lsqminnorm(AP_B,-A);
        E=lsqminnorm(AP_B,-C);
    end
    F=Y;
    X=E-P0;
elseif method==2
    X=-A*P0(npred+1:ndynamic,:);
    Y=[X-B(:,1:nspred),-B(:,nspred+1:ndynamic)];
    E=-C;
    F=-A;
    X_Y=B;
    if rank(B)<ndynamic
        X_0=X;
        X_YE=lsqminnorm(X_Y,E);
        X_YF=lsqminnorm(X_Y,F);
        X=-F*X_YE(npred+1:ndynamic,:);
        X_Y=X_Y+[X zeros(ndynamic,nfwrd)]-[zeros(ndynamic,npred) E*X_YF(1:nspred,:)];
        X=X_0+X;
        E=E*X_YE(1:nspred,:);
        F=F*X_YF(npred+1:ndynamic,:);
        j=j+1;
    end
end
X_0=X;

% if rank(AX_B)<ndynamic %If the starting matrix gives a singular coeffcient matrix
%     X=lsqminnorm(AX_B,-C);
%     AX_B=[A*X(npred+1:end,:) zeros(ndynamic,nfwrd)]+B;
%     M=AX_B*X+C;
%     j=1; % Counter for number of iterations
%     X_0=X;
% end

while diff>tol % As long as the iterations haven't converged...
    if method==1
        EIYX=eye(ndynamic)-[Y*X(npred+1:ndynamic,:) zeros(ndynamic,nfwrd)];
        EIYX=[E zeros(ndynamic, nfwrd)]/EIYX;
        FIXY=eye(ndynamic)-[zeros(ndynamic,npred) X*Y(1:nspred,:)];
        FIXY=[zeros(ndynamic, npred) F]/FIXY;
        %X_temp=EIYX*X;
        %Y_temp=FIXY*Y;
        X=X+FIXY*X*E(1:nspred,:);
        Y=Y+EIYX*Y*F(npred+1:ndynamic,:);
        E=EIYX*E;
        F=FIXY*F;
        rcond_coeff=min(rcond(EIYX'),rcond(FIXY'));%(A'\B')'
    elseif method==2
        X_YE=X_Y\E;
        X_YF=X_Y\F;
        X=-F*X_YE(npred+1:ndynamic,:);
        X_Y=X_Y+[X zeros(ndynamic,nfwrd)]-[zeros(ndynamic,npred) E*X_YF(1:nspred,:)];
        X=X_0+X;
        E=E*X_YE(1:nspred,:);
        F=F*X_YF(npred+1:ndynamic,:);
        rcond_coeff=rcond(X_Y);
    end
    diff=convergence_criterion([],X,X_0,[],A,B,C,convergence_metric);           % determine convergence criterion
    %diff=norm(X-X_0,1)/norm(X,1);
    j=j+1;                               % advance the counter
    X_0=X;
    if rcond_coeff<eps
        restart=restart+1;
        if restart<=max_restart && method==1
        EIYX=lsqminnorm([E zeros(ndynamic, nfwrd)]',EIYX')';
        FIXY=lsqminnorm([zeros(ndynamic, npred) F]',FIXY')';
        X=X+FIXY*X*E(1:nspred,:);
        Y=Y+EIYX*Y*F(npred+1:ndynamic,:);
        E=EIYX*E;
        F=FIXY*F;
        elseif restart<=max_restart && method==2
        X_YE=lsqminnorm(X_Y,E);
        X_YF=lsqminnorm(X_Y,F);
        X=-F*X_YE(npred+1:ndynamic,:);
        X_Y=X_Y+[X zeros(ndynamic,nfwrd)]-[zeros(ndynamic,npred) E*X_YF(1:nspred,:)];
        X=X_0+X;
        E=E*X_YE(1:nspred,:);
        F=F*X_YF(npred+1:ndynamic,:);
            j=j+1;
        else
            disp('Coefficient breakdown')
            j=1/0;
            break
        end
    end
    if j>max_it
        disp('Maximum iterations reached')
        break
    end
end

if method==1
    P=[real(P0+X) zeros(ndynamic,nfwrd)];
elseif method==2
    P=[real(-([A*P0(npred+1:ndynamic,:)+X+B(:,1:nspred) B(:,nspred+1:ndynamic)])\C) zeros(ndynamic,nfwrd)];
end

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




