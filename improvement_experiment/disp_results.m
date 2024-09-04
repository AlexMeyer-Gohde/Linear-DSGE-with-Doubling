load('First_run.mat')

results=HMS_Results;
table=[results(2,[1 6 4 5]) 1];
table=[table; [ results(3:end,1)./results(2,1) results(3:end,[ 6 4 5 7])]];
disp(compose('%1.2g',table))