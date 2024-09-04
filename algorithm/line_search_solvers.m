function [tj] = line_search_solvers_deflate(alpha_hk,beta_hk,gamma_hk,method)
% tj returns a scalar determining the size of the Newton step
%
% INPUTS
%   dPj, A, M   [matrix]     matrices dPj, A, M                                   
%
%   method      [string]     a string determining which method is used for
%                            calculating tj
%                                   
%
% OUTPUTS
%   tj          [scalar]     a scalar determining the size of the Newton
%                            step
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

%calculate elements of the merit function, following Higham & Kim, 2001
%(referenced in our paper)


% this method follows the body of the paper
if strcmp(method,'optimize')
    %method 1
    
    
    fun=@(x) (gamma_hk*x^4-beta_hk*x^3+(alpha_hk+beta_hk)*x^2-2*alpha_hk*x+alpha_hk);    %define function to be minimized
    options = optimset('Display','none');
    tj=fminbnd(fun,0,2,options);  % interval midpoint as starting value
    
% this method...
elseif strcmp(method,'direct')
    % method 2
    
    a=4*gamma_hk;
    b=-3*beta_hk;
    c=2*(alpha_hk+beta_hk);
    d=-2*alpha_hk;
    p=(3*a*c-b^2)/(3*a^2);
    q=(2*b^3-9*a*b*c+27*a^2*d)/(27*a^3);
    discr=-4*p^3-27*q^2;
    if discr>0+eps %three real roots
        %         if (a+b+c+d)*a>0 && a*b>0%No root or two roots in interval 0, inf % No root in interval 0, inf
        %             tj=2;
        %         else
        k=0:2;
        root=2*(-p/3)^(1/2)*cos(1/3*acos(3/2*q/p*(-3/p)^(1/2))-2*pi*k/3);
        root=root-b/(3*a);
        root=root(root>0&root<=2); root=[root,2];
        values=(gamma_hk*root.^4-beta_hk*root.^3+(alpha_hk+beta_hk)*root.^2-2*alpha_hk*root+alpha_hk);
        tj=root(values==min(values));
        %         end
    elseif discr<0-eps %one real root
        if p>0
            t=-2*(p/3)^(1/2)*sinh(1/3*asinh(3/2*q/p*(3/p)^(1/2)));
        else
            t=-2*abs(q)/q*(-p/3)^(1/2)*cosh(1/3*acosh(-3*abs(q)/(2*p)*(-3/p)^(1/2)));
        end
        %t=(-q/2+(q^2/4+p^3/27)^(1/2))^(1/3)+(-q/2-(q^2/4+p^3/27)^(1/2))^(1/3);
        %t=nthroot(-q/2+(q^2/4+p^3/27)^(1/2),3)+nthroot(-q/2-(q^2/4+p^3/27)^(1/2),3);
        tj=t-b/(3*a);
    else %discr=0 root =-b/(3*a)
        root=-b/(3*a);
        root=root(root>0&root<=2); root=[root,2];
        values=(gamma_hk*root.^4-beta_hk*root.^3+(alpha_hk+beta_hk)*root.^2-2*alpha_hk*root+alpha_hk);
        tj=root(values==min(values));
    end
    
else
    disp('Line search method not supported')
    tj=1;
end
%Pj=Pj+tj*dPj;                   % calculate solution matrix Pj
end

