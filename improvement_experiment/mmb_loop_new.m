%Description: ....
%....


clear 

run_time_reps=1;
%newton_options.dynare_reduced_sylvester=1;
newton_options.sylvester_method="dlyap_stripped";

%newton_options.sylvester_method="slicot";
newton_options.maximum_iterations=100;
addpath('C:\dynare\5.1\matlab')
addpath('..\algorithm\')


dynare Smets_Wouters_2007_std noclearall  nograph nostrict
%%%%% current problem: dynare_to_matrix_quadratic needs to be located in
%%%%% mmb-rep-folders (e.g. BRA_SAMBA08_rep)
HMS_Results(1,:)=[M_.nstatic, M_.nfwrd, M_.npred, M_.nboth, M_.nsfwrd, M_.nspred, M_.ndynamic];


[matrix_quadratic, jacobia_]=create_reduced_matrix_quadratic_from_dynare(M_,oo_);

%tic; [info, oo_, options_]  = stoch_simul(M_, options_, oo_, var_list_); toc    
total_time=[];  for jj=1:run_time_reps;tic;[dr,info] = dyn_first_order_solver(jacobia_,M_,oo_.dr,options_,0); total_time(jj)=toc;end;   HMS_Results(2,1) = mean(total_time(ceil(length(total_time)*.2):ceil(length(total_time)*.8)));   

 
 
ALPHA_ZS_dynare=[zeros(M_.endo_nbr,M_.nstatic) oo_.dr.ghx zeros(M_.endo_nbr,M_.nfwrd)];
X_dynare=ALPHA_ZS_dynare;
matrix_quadratic.X=ALPHA_ZS_dynare;
HMS_Results(2,6)=oo_.var(3,3);

 newton_options.initial=oo_.dr.ghx(M_.nstatic+1:end,:);
 bernoulli_options.initial=newton_options.initial;
try
[errors] = dsge_practical_forward_errors_matrix_quadratic(matrix_quadratic);
HMS_Results(2,3:5)=errors([4,7,8],1)';
catch
HMS_Results(2,3:5)=NaN(1,3);
end

try
 newton_options.algorithm='baseline';
 newton_options.M_=M_;
 newton_options.convergence_metric='reldiff';
 total_time=[]; for jj=1:run_time_reps;tic;[X,X_additional] = newton_matrix_quadratic(matrix_quadratic,newton_options); matrix_quadratic.X=X; [matrix_quadratic]=complete_reduced_matrix_quadratic(matrix_quadratic);total_time(jj)=toc;end;   
HMS_Results(3,1) =mean(total_time(ceil(length(total_time)*.2):ceil(length(total_time)*.8)));  
if max(max(isnan(matrix_quadratic.X)))==0; HMS_Results(3,2)=max(max(abs(X_dynare-matrix_quadratic.X)));end
HMS_Results(3,end)=X_additional;
oo.dr=oo_.dr;
oo.dr.ghu=[matrix_quadratic.B_full(:,1:M_.nstatic) [matrix_quadratic.A_static; matrix_quadratic.AA]*X(M_.npred+1:end,1:M_.nspred)+matrix_quadratic.B_full(:,M_.nstatic+1:end-M_.nfwrd) matrix_quadratic.B_full(:,end-M_.nfwrd+1:end)]; oo.dr.ghu=-oo.dr.ghu\matrix_quadratic.D; %oo.dr.ghu=-(matrix_quadratic.A_full*matrix_quadratic.X+matrix_quadratic.B_full)\matrix_quadratic.D;
oo.dr.ghx=matrix_quadratic.X(:,M_.nstatic+1:end-M_.nfwrd);
oo.dr.gx=oo.dr.ghx(end-M_.nfwrd+1:end,:);
oo= disp_th_moments(oo.dr, var_list, M_, options_, oo);
HMS_Results(3,6)=oo.var(3,3);
oo_newton=oo;
try
[errors] = dsge_practical_forward_errors_matrix_quadratic(matrix_quadratic);
HMS_Results(3,3:5)=errors([4,7,8],1)';
catch
HMS_Results(3,3:5)=NaN(1,3);
end
if X_additional==newton_options.maximum_iterations; 
HMS_Results(3,:)=NaN(1,7); 
end
catch
HMS_Results(3,:)=NaN(1,7);
end

clear bernoulli_options
bernoulli_options=[1 0 0 0 0 0 0 1 0 4 100 matrix_quadratic.ndynamic*eps 0];%matrix_quadratic.ndynamic*
bernoulli_options(13)=100;
% X_0_bernoulli=zeros(M_.ndynamic,M_.nspred);%X_0_1;
% [X,X_additional]=bernoulli_matrix_quadratic_careful(matrix_quadratic,X_0_bernoulli,bernoulli_options); 
% 
% matrix_quadratic.X=X; [matrix_quadratic]=complete_reduced_matrix_quadratic(matrix_quadratic);max(max(abs(X_dynare-matrix_quadratic.X)))
% maximum_iterations_total=bernoulli_options(11)+1;
% while X_additional==bernoulli_options(11)+1%&&max(max(abs(X_dynare-matrix_quadratic.X)))>1e-5
%     if maximum_iterations_total<1000000
%     bernoulli_options(11)=2*bernoulli_options(11);
%         [X,X_additional]=bernoulli_matrix_quadratic_careful(matrix_quadratic,X_0_bernoulli,bernoulli_options);
%         matrix_quadratic.X=X; [matrix_quadratic]=complete_reduced_matrix_quadratic(matrix_quadratic);max(max(abs(X_dynare-matrix_quadratic.X)))
%         %X_additional=X_additional+X_additional_plus;
%         maximum_iterations_total=maximum_iterations_total+bernoulli_options(11)
%         X_0_bernoulli=X(:,1:M_.nspred);
%     else
%         matrix_quadratic.X=X; [matrix_quadratic]=complete_reduced_matrix_quadratic(matrix_quadratic);
%         if max(max(abs(X_dynare-matrix_quadratic.X)))<1e-5
%             break
%         else
%             bernoulli_options(11)=1;
%             break
%         end
%     end
% end
% if bernoulli_options(11)~=1
%    bernoulli_options(11)=maximum_iterations_total;
% end
% 
% bernoulli_options(11)=100*bernoulli_options(11);

bernoulli_options(11)=50000;


try
 %   bernoulli_options.baseline=1;
    X_0_bernoulli=oo_.dr.ghx(M_.nstatic+1:end,:);%X_0_2;
total_time=[];  for jj=1:run_time_reps;tic;[X,X_additional]=bernoulli_matrix_quadratic_careful(matrix_quadratic,X_0_bernoulli,bernoulli_options); matrix_quadratic.X=X; [matrix_quadratic]=complete_reduced_matrix_quadratic(matrix_quadratic); total_time(jj)=toc;end; 
HMS_Results(4,1) = mean(total_time(ceil(length(total_time)*.2):ceil(length(total_time)*.8)));  
if max(max(isnan(matrix_quadratic.X)))==0; HMS_Results(4,2)=max(max(abs(X_dynare-matrix_quadratic.X)));end
HMS_Results(4,end)=X_additional;
X_0_2=X(:,1:M_.nspred);
oo.dr=oo_.dr;
oo.dr.ghu=[matrix_quadratic.B_full(:,1:M_.nstatic) [matrix_quadratic.A_static; matrix_quadratic.AA]*X(M_.npred+1:end,1:M_.nspred)+matrix_quadratic.B_full(:,M_.nstatic+1:end-M_.nfwrd) matrix_quadratic.B_full(:,end-M_.nfwrd+1:end)]; oo.dr.ghu=-oo.dr.ghu\matrix_quadratic.D; %oo.dr.ghu=-(matrix_quadratic.A_full*matrix_quadratic.X+matrix_quadratic.B_full)\matrix_quadratic.D;
oo.dr.ghx=matrix_quadratic.X(:,M_.nstatic+1:end-M_.nfwrd);
oo.dr.gx=oo.dr.ghx(end-M_.nfwrd+1:end,:);
oo= disp_th_moments(oo.dr, var_list, M_, options_, oo);
HMS_Results(4,6)=oo.var(3,3);
try
[errors] = dsge_practical_forward_errors_matrix_quadratic(matrix_quadratic);
HMS_Results(4,3:5)=errors([4,7,8],1)';
catch
HMS_Results(4,3:5)=NaN(1,3);
end
% if X_additional==newton_options.maximum_iterations; 
% HMS_Results(4,:)=NaN(1,7); 
% end
catch
HMS_Results(4,:)=NaN(1,7);
end
oo_bernoulli=oo;

%bernoulli_options.maximum_iterations=1;


try
 
doubling_options.method='SF1';
doubling_options.P0=oo_.dr.ghx(M_.nstatic+1:end,:);
total_time=[];  for jj=1:run_time_reps;tic;[X,X_additional]=doubling_matrix_quadratic(matrix_quadratic,doubling_options); matrix_quadratic.X=X; [matrix_quadratic]=complete_reduced_matrix_quadratic(matrix_quadratic); total_time(jj)=toc;end; 
HMS_Results(5,1) = mean(total_time(ceil(length(total_time)*.2):ceil(length(total_time)*.8))); 
if max(max(isnan(matrix_quadratic.X)))==0; HMS_Results(5,2)=max(max(abs(X_dynare-matrix_quadratic.X)));end
HMS_Results(5,end)=X_additional;
X_0_3=X(:,1:M_.nspred);
oo.dr=oo_.dr;
oo.dr.ghu=[matrix_quadratic.B_full(:,1:M_.nstatic) [matrix_quadratic.A_static; matrix_quadratic.AA]*X(M_.npred+1:end,1:M_.nspred)+matrix_quadratic.B_full(:,M_.nstatic+1:end-M_.nfwrd) matrix_quadratic.B_full(:,end-M_.nfwrd+1:end)]; oo.dr.ghu=-oo.dr.ghu\matrix_quadratic.D; %oo.dr.ghu=-(matrix_quadratic.A_full*matrix_quadratic.X+matrix_quadratic.B_full)\matrix_quadratic.D;
oo.dr.ghx=matrix_quadratic.X(:,M_.nstatic+1:end-M_.nfwrd);
oo.dr.gx=oo.dr.ghx(end-M_.nfwrd+1:end,:);
oo= disp_th_moments(oo.dr, var_list, M_, options_, oo);
HMS_Results(5,6)=oo.var(3,3);
try
[errors] = dsge_practical_forward_errors_matrix_quadratic(matrix_quadratic);
HMS_Results(5,3:5)=errors([4,7,8],1)';
catch
HMS_Results(5,3:5)=NaN(1,3);
end
% if X_additional==newton_options.maximum_iterations; 
% HMS_Results(5,:)=NaN(1,7); 
% end
catch
HMS_Results(5,:)=NaN(1,7);
end

try

doubling_options.method='SF2';
doubling_options.P0=oo_.dr.ghx(M_.nstatic+1:end,:);
total_time=[];  for jj=1:run_time_reps;tic;[X,X_additional]=doubling_matrix_quadratic(matrix_quadratic,doubling_options); matrix_quadratic.X=X; [matrix_quadratic]=complete_reduced_matrix_quadratic(matrix_quadratic); total_time(jj)=toc;end; 
HMS_Results(6,1) = mean(total_time(ceil(length(total_time)*.2):ceil(length(total_time)*.8)));
if max(max(isnan(matrix_quadratic.X)))==0; HMS_Results(6,2)=max(max(abs(X_dynare-matrix_quadratic.X)));end
HMS_Results(6,end)=X_additional;
X_0_4=X(:,1:M_.nspred);
oo.dr=oo_.dr;
oo.dr.ghu=[matrix_quadratic.B_full(:,1:M_.nstatic) [matrix_quadratic.A_static; matrix_quadratic.AA]*X(M_.npred+1:end,1:M_.nspred)+matrix_quadratic.B_full(:,M_.nstatic+1:end-M_.nfwrd) matrix_quadratic.B_full(:,end-M_.nfwrd+1:end)]; oo.dr.ghu=-oo.dr.ghu\matrix_quadratic.D; %oo.dr.ghu=-(matrix_quadratic.A_full*matrix_quadratic.X+matrix_quadratic.B_full)\matrix_quadratic.D;
oo.dr.ghx=matrix_quadratic.X(:,M_.nstatic+1:end-M_.nfwrd);
oo.dr.gx=oo.dr.ghx(end-M_.nfwrd+1:end,:);
oo= disp_th_moments(oo.dr, var_list, M_, options_, oo);
HMS_Results(6,6)=oo.var(3,3);
try
[errors] = dsge_practical_forward_errors_matrix_quadratic(matrix_quadratic);
HMS_Results(6,3:5)=errors([4,7,8],1)';
catch
HMS_Results(6,3:5)=NaN(1,3);
end
% if X_additional==newton_options.maximum_iter12,ations; 
% HMS_Results(6,:)=NaN(1,7); 
% end
catch
HMS_Results(6,:)=NaN(1,7);
end

%HMS_Results(3:end,6)=(HMS_Results(3:end,6)-HMS_Results(2,6))/HMS_Results(2,6);


%run different Newton methods
%newton_solvent_2(matrix_quadratic);
%newton_solvent_3(matrix_quadratic);
 save First_run


