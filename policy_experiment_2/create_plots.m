% clear all
% load('Policy_Run.mat')
% 
% combined_results_old=combined_results([1,2,7,8,3,4,5,6],:,:);
% 
% load('Policy_Run_fill.mat')
% 
% combined_results_old(7:8,:,:)=combined_results(5:6,:,:);
% combined_results=combined_results_old;

% clear all
% load('Policy_Run_2.mat')
% combined_results(3:8,:,:)=combined_results([7:8,3:6],:,:);

% clear all
% load('Policy_Run_2.mat')
% combined_results_old=combined_results([1,2,7,8,3,4,5,6],:,:);
% load('Policy_Run_2_fill.mat')
% combined_results_old([5 7 8],:,:)=combined_results([3 5 6],:,:);
% load('Policy_Run_2_fill.mat')
% combined_results_old([5],:,:)=combined_results([3],:,:);
% combined_results=combined_results_old;

clear all
load('Policy_Run_5.mat')
combined_results(3:8,:,:)=combined_results([7:8,3:6],:,:);

grid_size=STEPS;
grid_size=10.^grid_size;
highlight=[1:7];
newcolors=[0 0 0;
    0 0 1
1 0 0
0.929000000000000	0.694000000000000	0.125000000000000
0.494000000000000	0.184000000000000	0.556000000000000
0.466000000000000	0.674000000000000	0.188000000000000
0.301000000000000	0.745000000000000	0.933000000000000
0.933000000000000	0.000000000000000	0.933000000000000];
newcolors=newcolors([1 2 3 4 5 6 7],:);
lineorder={'-','-','-','-','-','-','-',':','-','--',':','-','--',};
colororder(newcolors);
figure
set(gcf,'DefaultAxesColorOrder',newcolors,'DefaultAxesLineStyleOrder',lineorder)
ax = gca; 
hold on
plot(log10(grid_size),log10(squeeze(combined_results(2,1,:))),'LineWidth',2)
  ax.LineStyleOrderIndex = ax.ColorOrderIndex;
for j=2:7
   if sum(highlight==j)
plot(log10(grid_size),log10(squeeze(combined_results(j+1,1,:))),'LineWidth',2)
   else
plot(log10(grid_size),log10(squeeze(combined_results(j+1,1,:))),'LineWidth',1)
   end
     ax.LineStyleOrderIndex = ax.ColorOrderIndex;
end
% hold off
% yl = ylim
% ylim=([yl(1)-1 yl(2)]);
% ylim
legend('QZ','Dynare (CR)','Dynare (LR)','Newton','Bernoulli','SF1','SF2','AutoUpdate','off','location','south','NumColumns',4)
ylabel('Computation Time per Grid Point, Seconds, Log10')
xlabel('Grid Point Distance, -Log10')
'Dynare (CR)','Dynare (LR)','Newton','Bernoulli'
%colororder(newcolors);
%C = colororder;
%newcolors=C(2:end,:);

newcolors=newcolors(2:end,:);
lineorder=lineorder(2:end);
% figure
% hold on
% for j=2:10
% scatter(log10(squeeze(combined_results(j+1,4,:))./squeeze(combined_results(2,4,:))),log10(squeeze(combined_results(j+1,1,:))./squeeze(combined_results(2,1,:))),'filled')
% end
% legend('Baseline','Modified','Samanskii','LS', 'Occ. Line','Occ. Line Sam','AutoUpdate','off','location','southeast')
% ylabel('Computation Time, Relative to Dynare, Log10')
% xlabel('Forward Error Bound 1, Relative to Dynare, Log10')
% xlim([-1.2 .2])
% %ylim([-2 2])
% %set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin')
% xline(0);
% yline(0);
% hold off
% 
% figure
% hold on
% for j=2:10
% scatter(log10(squeeze(combined_results(j+1,5,:))./squeeze(combined_results(2,5,:))),log10(squeeze(combined_results(j+1,1,:))./squeeze(combined_results(2,1,:))),'filled')
% end
% legend('Baseline','Modified','Samanskii','LS', 'Occ. Line','Occ. Line Sam','AutoUpdate','off','location','southeast')
% ylabel('Computation Time, Relative to Dynare, Log10')
% xlabel('Forward Error Bound 2, Relative to Dynare, Log10')
% xlim([-1.2 .2])
% %ylim([-2 2])
% %set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin')
% xline(0);
% yline(0);
% hold off

figure
set(gcf,'DefaultAxesColorOrder',newcolors,'DefaultAxesLineStyleOrder',lineorder)
ax = gca; 
hold on
for j=2:7
   if sum(highlight==j)
plot(log10(grid_size),log10(squeeze(combined_results(j+1,4,:))./squeeze(combined_results(2,4,:))),'LineWidth',2)
   else
plot(log10(grid_size),log10(squeeze(combined_results(j+1,4,:))./squeeze(combined_results(2,4,:))),'LineWidth',1)
   end
     ax.LineStyleOrderIndex = ax.ColorOrderIndex;
end
legend('Dynare (CR)','Dynare (LR)','Newton','Bernoulli','SF1','SF2','AutoUpdate','off','location','southeast')
ylabel('Forward Error Bound 1, Relative to Dynare, Log10')
xlabel('Grid Point Distance, -Log10')
colororder(newcolors);

figure
set(gcf,'DefaultAxesColorOrder',newcolors,'DefaultAxesLineStyleOrder',lineorder)
ax = gca; 
hold on
for j=2:7
   if sum(highlight==j)
plot(log10(grid_size),log10(squeeze(combined_results(j+1,5,:))./squeeze(combined_results(2,5,:))),'LineWidth',2)
   else
plot(log10(grid_size),log10(squeeze(combined_results(j+1,5,:))./squeeze(combined_results(2,5,:))),'LineWidth',1)
   end
     ax.LineStyleOrderIndex = ax.ColorOrderIndex;
end
legend('Dynare (CR)','Dynare (LR)','Newton','Bernoulli','SF1','SF2','AutoUpdate','off','location','southeast')
ylabel('Forward Error Bound 2, Relative to Dynare, Log10')
xlabel('Grid Point Distance, -Log10')
colororder(newcolors);

newcolors=[0 0 0;newcolors];
lineorder={'-',lineorder{:}};
figure
set(gcf,'DefaultAxesColorOrder',newcolors,'DefaultAxesLineStyleOrder',lineorder)
ax = gca; 
hold on
plot(log10(grid_size),log10(squeeze(combined_results(2,4,:))),'LineWidth',2)
for j=2:7
   if sum(highlight==j)
plot(log10(grid_size),log10(squeeze(combined_results(j+1,4,:))),'LineWidth',2)
   else
plot(log10(grid_size),log10(squeeze(combined_results(j+1,4,:))),'LineWidth',1)
   end
     ax.LineStyleOrderIndex = ax.ColorOrderIndex;
end
legend('QZ','Dynare (CR)','Dynare (LR)','Newton','Bernoulli','SF1','SF2','AutoUpdate','off','location','southeast')
ylabel('Forward Error Bound 1,  Log10')
xlabel('Grid Point Distance, -Log10')
colororder(newcolors);

figure
set(gcf,'DefaultAxesColorOrder',newcolors,'DefaultAxesLineStyleOrder',lineorder)
ax = gca; 
hold on
plot(log10(grid_size),log10(squeeze(combined_results(2,5,:))),'LineWidth',2)
for j=2:7
    if sum(highlight==j)
plot(log10(grid_size),log10(squeeze(combined_results(j+1,5,:))),'LineWidth',2)
    else
plot(log10(grid_size),log10(squeeze(combined_results(j+1,5,:))),'LineWidth',1)
    end
     ax.LineStyleOrderIndex = ax.ColorOrderIndex;
end
legend('QZ','Dynare (CR)','Dynare (LR)','Newton','Bernoulli','SF1','SF2','AutoUpdate','off','location','southeast')
ylabel('Forward Error Bound 2,  Log10')
xlabel('Grid Point Distance, -Log10')
colororder(newcolors);