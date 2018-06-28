clear all;clc;
ratedata = csvread('FBFitted.csv',1,0);
[m,n] = size(ratedata);
yielddata = log(1.+ ratedata(:,2:n)./100);
yielddata = [ratedata(:,1),yielddata];
temmat = ones(31,n);
timeseries = [1:31]';
for i = 1:31
    for j = 1:m
        if j/12 == i
          temmat(i,:) = yielddata(j,:);
        end
    end
end
annualdata = [timeseries,temmat(:,[6,10,12:19])];
cons = ones(30,1);
%% Question one
display('Question one');
Beta1 = zeros(1,2);
for n = 2:3
    Expectdi_y = annualdata([2:31],n) - annualdata([1:30],n+1);
    di_x = (annualdata([1:30],n+1) - annualdata([1:30],2))/(n-1);
    [A1,BINT1,R1,RINT1,STATS1] = regress(Expectdi_y,[di_x]);
    fprintf('For %1.0f year bonds, Beta is %3.4f.\n',n,A1);
end
%% Question Two
display('Question two');
n = 2;
Rm = n*annualdata([1:30],n+1)-(n-1)*annualdata([2:31],n);
Sm = (n-1)*(annualdata([1:30],n+1)-annualdata([2:31],n))-annualdata([1:30],n+1)-annualdata([1:30],2);
Mdl = varm(2,1);
EstMdl = estimate(Mdl,[Rm Sm]);
summarize(EstMdl);
AR = cell2mat(EstMdl.AR);
shock = [Rm(2:30) Sm(2:30)]'-EstMdl.Constant-AR*[Rm(1:29) Sm(1:29)]';
out_Rm = decomp(Rm(2:30),shock(1,:)',shock(2,:)')
out_Sm = decomp(Sm(2:30),shock(1,:)',shock(2,:)')
%% Question Three
display('Question three');
Beta3 = zeros(1,9);
p_value3 = zeros(1,9);
for n = 2:10
    Expectdi = n*annualdata([1:30],n+1)-(n-1)*annualdata([2:31],n)-annualdata([1:30],2);
    di = n*annualdata([1:30],n+1)-(n-1)*annualdata([1:30],n)-annualdata([1:30],2);
    stats = regstats(Expectdi,di);
    Beta3(n-1) = stats.beta(2);
    p_value3(n-1) = stats.tstat.pval(2);
end
for n = 2:10
    if p_value3(n-1) <= 0.05
    fprintf('For %1.0f year bonds, Beta is %5.4f and p_value is %5.4f, and the expected hypothesis is wrong. \n',...
        n,Beta3(n-1),p_value3(n-1));
    else
        fprintf('For %1.0f year bonds, Beta is %5.4f and p_value is %5.4f, and the expected hypothesis is ture. \n',...
        n,Beta3(n-1),p_value3(n-1));
    end
end