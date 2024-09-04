
load('First_Run_HMS_improvement.mat')
HMS_Results_temp=HMS_Results;
load('First_Run_HMS_improvement_fill.mat')
HMS_Results_temp(5:6,:,:)=HMS_Results(5:6,:,:);
HMS_Results=HMS_Results_temp;


model_variable_number=sum(squeeze(HMS_Results(1,1:4,:)));

times=squeeze(HMS_Results(2:end,1,:));
iterations=squeeze(HMS_Results(2:end,7,:)); iterations(1,:)=1;
rel_times=times./times(1,:);
diffs=squeeze(HMS_Results(2:end,2,:));
diffs_select=0<diffs&diffs<1e-5;
diffs_select(1,:)=1;

fe_1=squeeze(HMS_Results(2:end,4,:));
rel_fe_1=fe_1./fe_1(1,:);
fe_2=squeeze(HMS_Results(2:end,5,:));
rel_fe_2=fe_2./fe_2(1,:);



results(1:5,1)=sum(diffs_select')';
for j=1:5
results(j,2)=median(rel_times(j,diffs_select(j,:)),'omitnan');
results(j,3)=min(rel_times(j,diffs_select(j,:)));
results(j,4)=max(rel_times(j,diffs_select(j,:)));
results(j,5)=median(rel_fe_1(j,diffs_select(j,:)),'omitnan');
results(j,6)=min(rel_fe_1(j,diffs_select(j,:)));
results(j,7)=max(rel_fe_1(j,diffs_select(j,:)));
results(j,8)=median(rel_fe_2(j,diffs_select(j,:)),'omitnan');
results(j,9)=min(rel_fe_2(j,diffs_select(j,:)));
results(j,10)=max(rel_fe_2(j,diffs_select(j,:)));
results(j,11)=median(iterations(j,diffs_select(j,:)),'omitnan');
end

disp(compose('%1.2g',results))
compose('%1.2g',results)


highlight=[1:7];
newcolors=[0 0 0;
%    0 0 1
%1 0 0
0.929000000000000	0.694000000000000	0.125000000000000
0.494000000000000	0.184000000000000	0.556000000000000
0.466000000000000	0.674000000000000	0.188000000000000
0.301000000000000	0.745000000000000	0.933000000000000
0.933000000000000	0.000000000000000	0.933000000000000];
newcolors=newcolors([1 2 3 4 5 ],:);
lineorder={'-','-','-','-','-','-','-',':','-','--',':','-','--',};
colororder(newcolors);

figure
set(gcf,'DefaultAxesColorOrder',newcolors,'DefaultAxesLineStyleOrder',lineorder)
ax = gca; 
hold on
for j=2:5
       if sum(highlight==j)
           scatter(log10(rel_fe_1(j,diffs_select(j,:))),log10(rel_times(j,diffs_select(j,:))),'LineWidth',1)
           flag=0;
       else
           if flag==0
           scatter(log10(rel_fe_1(j,diffs_select(j,:))),log10(rel_times(j,diffs_select(j,:))),'x','LineWidth',1)
           flag=1;
           else
           scatter(log10(rel_fe_1(j,diffs_select(j,:))),log10(rel_times(j,diffs_select(j,:))),'+','LineWidth',1)
           end
       end
%scatter(log10(model_variable_number(diffs_select(j,:))),log10(rel_fe_2(j,diffs_select(j,:))),'filled')
ax.LineStyleOrderIndex = ax.ColorOrderIndex;
end
legend('Newton','Bernoulli','SF1','SF2','AutoUpdate','off','location','southeast')
ylabel('Computation Time, Relative to Dynare, Log10')
xlabel('Forward Error Bound 1, Relative to Dynare, Log10')
%xlim([-2 2])
%ylim([-2 2])
%set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin')
xline(0);
yline(0);
colororder(newcolors);
hold off

legends={'Newton','Bernoulli','SF1','SF2'};

for j=2:5
    hold on
    figure
scatter(log10(rel_fe_1(j,diffs_select(j,:))),log10(rel_times(j,diffs_select(j,:))),'filled')
legend(legends{j-1},'AutoUpdate','off','location','northeast')
%legend('Baseline','Modified','Samanskii','Line Search', 'Occ. Line','Occ. Line Search Samanskii','AutoUpdate','off','location','southeast')
ylabel('Computation Time, Relative to Dynare, Log10')
xlabel('Forward Error Bound 1, Relative to Dynare, Log10')
xlim([-4 8])
ylim([-1 3])
h=lsline;
h.Color = 'b';
h.LineStyle = '--';
%set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin')
xline(0);
yline(0);
hold off
end

figure
set(gcf,'DefaultAxesColorOrder',newcolors,'DefaultAxesLineStyleOrder',lineorder)
ax = gca; 
hold on
for j=2:5
       if sum(highlight==j)
           scatter(log10(rel_fe_2(j,diffs_select(j,:))),log10(rel_times(j,diffs_select(j,:))),'LineWidth',1)
           flag=0;
       else
           if flag==0
           scatter(log10(rel_fe_2(j,diffs_select(j,:))),log10(rel_times(j,diffs_select(j,:))),'x','LineWidth',1)
           flag=1;
           else
           scatter(log10(rel_fe_2(j,diffs_select(j,:))),log10(rel_times(j,diffs_select(j,:))),'+','LineWidth',1)
           end
       end
%scatter(log10(model_variable_number(diffs_select(j,:))),log10(rel_fe_2(j,diffs_select(j,:))),'filled')
ax.LineStyleOrderIndex = ax.ColorOrderIndex;
end
legend('Newton','Bernoulli','SF1','SF2','AutoUpdate','off','location','southeast')
ylabel('Computation Time, Relative to Dynare, Log10')
xlabel('Forward Error Bound 2, Relative to Dynare, Log10')
%xlim([-2 2])
%ylim([-2 2])
%set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin')
xline(0);
yline(0);
colororder(newcolors);
hold off


for j=2:5
    hold on
    figure
scatter(log10(rel_fe_2(j,diffs_select(j,:))),log10(rel_times(j,diffs_select(j,:))),'filled');
legend(legends{j-1},'AutoUpdate','off','location','northeast')

%legend('Baseline','Modified','Samanskii','Line Search', 'Occ. Line Search','Occ. Line Sam','AutoUpdate','off','location','southeast')
ylabel('Computation Time, Relative to Dynare, Log10')
xlabel('Forward Error Bound 2, Relative to Dynare, Log10')
xlim([-6 6])
ylim([-1 3])
h=lsline;
h.Color = 'b';
h.LineStyle = '--';
%set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin')
xline(0);
yline(0);
hold off
end




figure
set(gcf,'DefaultAxesColorOrder',newcolors,'DefaultAxesLineStyleOrder',lineorder)
ax = gca; 
hold on
for j=1:5
    [y,x]=ksdensity(log10(fe_1(j,diffs_select(j,:)))); 
    if sum(highlight==j)
        plot(x,y,'LineWidth',2)
    else
        plot(x,y,'LineWidth',1)
    end
     ax.LineStyleOrderIndex = ax.ColorOrderIndex;
end
legend('QZ','Newton','Bernoulli','SF1','SF2','AutoUpdate','off','location','southeast')
ylabel('Density')
xlabel('Forward Error Bound 1,  Log10')
xlim([-18 -10])
%ylim([-2 2])
%set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin')
xl =xline(log10(eps),'-.','Machine Precision','LineWidth',2);
%xl.LabelVerticalAlignment = 'top';
xl.LabelHorizontalAlignment = 'left';
legend('show');
hold off

figure
set(gcf,'DefaultAxesColorOrder',newcolors,'DefaultAxesLineStyleOrder',lineorder)
ax = gca; 
hold on
for j=1:5
    [y,x]=ksdensity(log10(fe_2(j,diffs_select(j,:)))); 
    if sum(highlight==j)
        plot(x,y,'LineWidth',2)
    else
        plot(x,y,'LineWidth',1)
    end
     ax.LineStyleOrderIndex = ax.ColorOrderIndex;
end
legend('QZ','Newton','Bernoulli','SF1','SF2','AutoUpdate','off','location','southeast')
ylabel('Density')
xlabel('Forward Error Bound 2,  Log10')
xlim([-18 0])
%ylim([-2 2])
%set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin')
xl =xline(log10(eps),'-.','Machine Precision','LineWidth',2);
%xl.LabelVerticalAlignment = 'top';
xl.LabelHorizontalAlignment = 'left';
legend('show');
hold off


newcolors=newcolors(2:end,:);
lineorder=lineorder(2:end);

figure
set(gcf,'DefaultAxesColorOrder',newcolors,'DefaultAxesLineStyleOrder',lineorder)
ax = gca; 
hold on
for j=2:5
    [y,x]=ksdensity(log10(rel_fe_1(j,diffs_select(j,:)))); 
    if sum(highlight==j)
        plot(x,y,'LineWidth',2)
    else
        plot(x,y,'LineWidth',1)
    end
     ax.LineStyleOrderIndex = ax.ColorOrderIndex;
end
legend('Newton','Bernoulli','SF1','SF2','AutoUpdate','off','location','southeast')
ylabel('Density')
xlabel('Forward Error Bound 1, Relative to Dynare, Log10')
xlim([-3 3])
%ylim([-2 2])
%set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin')
xl =xline(log10(1),'-.','LineWidth',2);
%xl.LabelVerticalAlignment = 'top';
%xl.LabelHorizontalAlignment = 'left';
%legend('show');
hold off

figure
set(gcf,'DefaultAxesColorOrder',newcolors,'DefaultAxesLineStyleOrder',lineorder)
ax = gca; 
hold on
for j=2:5
    [y,x]=ksdensity(log10(rel_fe_2(j,diffs_select(j,:)))); 
    if sum(highlight==j)
        plot(x,y,'LineWidth',2)
    else
        plot(x,y,'LineWidth',1)
    end
     ax.LineStyleOrderIndex = ax.ColorOrderIndex;
end
legend('Newton','Bernoulli','SF1','SF2','AutoUpdate','off','location','southeast')
ylabel('Density')
xlabel('Forward Error Bound 2, Relative to Dynare, Log10')
xlim([-3 3])
%ylim([-2 2])
%set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin')
xl =xline(log10(1),'-.','LineWidth',2);
%xl.LabelVerticalAlignment = 'top';
%xl.LabelHorizontalAlignment = 'left';
%legend('show');
hold off

figure
set(gcf,'DefaultAxesColorOrder',newcolors,'DefaultAxesLineStyleOrder',lineorder)
ax = gca; 
hold on
for j=2:5
           if sum(highlight==j)
scatter(log10(rel_fe_1(j,diffs_select(j,:))),log10(rel_times(j,diffs_select(j,:))),'LineWidth',1)
           flag=0;
           else
           if flag==0
               scatter(log10(rel_fe_1(j,diffs_select(j,:))),log10(rel_times(j,diffs_select(j,:))),'x','LineWidth',1)
                          flag=1;
           else
               scatter(log10(rel_fe_1(j,diffs_select(j,:))),log10(rel_times(j,diffs_select(j,:))),'+','LineWidth',1)
           end
           end
     ax.LineStyleOrderIndex = ax.ColorOrderIndex;
end
legend('Newton','Bernoulli','SF1','SF2','AutoUpdate','off','location','southeast')
ylabel('Computation Time, Relative to Dynare, Log10')
xlabel('Forward Error Bound 1, Relative to Dynare, Log10')
%xlim([-2 2])
%ylim([-2 2])
%set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin')
xline(0);
yline(0);
hold off

figure
set(gcf,'DefaultAxesColorOrder',newcolors,'DefaultAxesLineStyleOrder',lineorder)
ax = gca; 
hold on
for j=2:5
           if sum(highlight==j)
scatter(log10(rel_fe_2(j,diffs_select(j,:))),log10(rel_times(j,diffs_select(j,:))),'LineWidth',1)
           flag=0;
           else
           if flag==0
               scatter(log10(rel_fe_2(j,diffs_select(j,:))),log10(rel_times(j,diffs_select(j,:))),'x','LineWidth',1)
                          flag=1;
           else
               scatter(log10(rel_fe_2(j,diffs_select(j,:))),log10(rel_times(j,diffs_select(j,:))),'+','LineWidth',1)
           end
           end
     ax.LineStyleOrderIndex = ax.ColorOrderIndex;
end
legend('Newton','Bernoulli','SF1','SF2','AutoUpdate','off','location','southeast')
ylabel('Computation Time, Relative to Dynare, Log10')
xlabel('Forward Error Bound 2, Relative to Dynare, Log10')
%xlim([-2 2])
%ylim([-2 2])
%set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin')
xline(0);
yline(0);
hold off

%Plot computation time relative to dynare against model size (paper:
%figure x)
figure
set(gcf,'DefaultAxesColorOrder',newcolors,'DefaultAxesLineStyleOrder',lineorder)
ax = gca; 
hold on
% for j=2:5
% scatter(log10(model_variable_number(diffs_select(j,:))),log10(rel_times(j,diffs_select(j,:))),'filled')
% end
for j=2:5
       if sum(highlight==j)
           scatter(log10(model_variable_number(diffs_select(j,:))),log10(rel_times(j,diffs_select(j,:))),'LineWidth',1)
           flag=0;
       else
           if flag==0
           scatter(log10(model_variable_number(diffs_select(j,:))),log10(rel_times(j,diffs_select(j,:))),'x','LineWidth',1)
           flag=1;
           else
           scatter(log10(model_variable_number(diffs_select(j,:))),log10(rel_times(j,diffs_select(j,:))),'+','LineWidth',1)
           end
       end
%scatter(log10(model_variable_number(diffs_select(j,:))),log10(rel_fe_2(j,diffs_select(j,:))),'filled')
ax.LineStyleOrderIndex = ax.ColorOrderIndex;
end
legend('Newton','Bernoulli','SF1','SF2','AutoUpdate','off','location','southeast')
ylabel('Computation Time, Relative to Dynare, Log10')
xlabel('Model Size, Log10')
xline(0);
yline(0);
hold off

%Plot forward error bounds relative to dynare against model size (paper:
%figure x)
figure
set(gcf,'DefaultAxesColorOrder',newcolors,'DefaultAxesLineStyleOrder',lineorder)
ax = gca; 
hold on
%for j=2:5
%scatter(log10(model_variable_number(diffs_select(j,:))),log10(rel_fe_1(j,diffs_select(j,:))),'filled')
%end
for j=2:5
       if sum(highlight==j)
           scatter(log10(model_variable_number(diffs_select(j,:))),log10(rel_fe_1(j,diffs_select(j,:))),'LineWidth',1)
           flag=0;
       else
           if flag==0
           scatter(log10(model_variable_number(diffs_select(j,:))),log10(rel_fe_1(j,diffs_select(j,:))),'x','LineWidth',1)
           flag=1;
           else
           scatter(log10(model_variable_number(diffs_select(j,:))),log10(rel_fe_1(j,diffs_select(j,:))),'+','LineWidth',1)
           end
       end
%scatter(log10(model_variable_number(diffs_select(j,:))),log10(rel_fe_2(j,diffs_select(j,:))),'filled')
ax.LineStyleOrderIndex = ax.ColorOrderIndex;
end
legend('Newton','Bernoulli','SF1','SF2','AutoUpdate','off','location','southeast')
ylabel('Forward Error Bound 1, Relative to Dynare, Log10')
xlabel('Model Size, Log10')
xline(0);
yline(0);
ylim([-6 8])
colororder(newcolors);
hold off



figure
set(gcf,'DefaultAxesColorOrder',newcolors,'DefaultAxesLineStyleOrder',lineorder)
ax = gca; 
hold on
for j=2:5
       if sum(highlight==j)
           scatter(log10(model_variable_number(diffs_select(j,:))),log10(rel_fe_2(j,diffs_select(j,:))),'LineWidth',1)
           flag=0;
       else
           if flag==0
           scatter(log10(model_variable_number(diffs_select(j,:))),log10(rel_fe_2(j,diffs_select(j,:))),'x','LineWidth',1)
           flag=1;
           else
           scatter(log10(model_variable_number(diffs_select(j,:))),log10(rel_fe_2(j,diffs_select(j,:))),'+','LineWidth',1)
           end
       end
%scatter(log10(model_variable_number(diffs_select(j,:))),log10(rel_fe_2(j,diffs_select(j,:))),'filled')
ax.LineStyleOrderIndex = ax.ColorOrderIndex;
end
legend('Newton','Bernoulli','SF1','SF2','AutoUpdate','off','location','southeast')
ylabel('Forward Error Bound 2, Relative to Dynare, Log10')
xlabel('Model Size, Log10')
colororder(newcolors);
xline(0);
yline(0);
ylim([-6 8])
hold off





for j=2:5
    hold on
    figure
scatter(log10(model_variable_number(diffs_select(j,:))),log10(rel_fe_1(j,diffs_select(j,:))),'filled');
legend(legends{j-1},'AutoUpdate','off','location','northeast')

%legend('Baseline','Modified','Samanskii','Line Search', 'Occ. Line Search','Occ. Line Sam','AutoUpdate','off','location','southeast')
ylabel('Forward Error Bound 1, Relative to Dynare, Log10')
xlabel('Model Size, Log10')
xlim([0 3.5])
%ylim([-6 8])
h=lsline;
h.Color = 'b';
h.LineStyle = '--';
%set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin')
xline(0);
yline(0);
hold off
end

disp(compose('%1.2g',results))
