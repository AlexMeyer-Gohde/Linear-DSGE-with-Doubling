function [diff] = convergence_criterion(M,P,P0,G,A,B,C,method)
% diff returns the convergence criterion
%
% INPUTS
%   M, Pj, P0, G,
%   A, B, C      [matrices]     matrices M,Pj,P0,G,A,B,C                               
%
%   method       [string]       A string determining which method is used 
%                               for calculating the convergence criterion   
%
% OUTPUTS
%   diff         [scalar]       The convergence criterion
%
%
% Copyright (C) 2022 Alexander Meyer-Gohde & Johanna Saecker
%
% This is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

% the default convergence criterion follows ...
if strcmp(method,'residual')
    diff=norm(M,'fro');%=norm((A*P+B)*P+C,'fro');
elseif strcmp(method,'diff')
    diff=norm(P-P0,1);
elseif strcmp(method,'reldiff')
    diff=norm(P-P0,1)/norm(P,1);
% elseif strcmp(method,'reldiffX')
%     diff=norm(P-P0,1)/norm(X,1);
elseif strcmp(method,'fe1')
    R_P=M;
    V=kron(eye(size(P,2)),G)+kron(P(1:size(C,2),:)',[zeros(size(B,2),size(B,2)-size(A,2)) A]);
    temp_P=lsqminnorm(V,R_P(:));
    diff=norm(temp_P)/norm(P,'fro');
elseif strcmp(method,'fe1_sparse')
% R_P=M;
% V=kron(speye(size(B,1)),G)+kron(P',A);
%     temp_P=lsqminnorm(V,R_P(:));
% diff=norm(temp_P)/P_F;
else
    disp('Unknown Convergence Criterion');
end
end

