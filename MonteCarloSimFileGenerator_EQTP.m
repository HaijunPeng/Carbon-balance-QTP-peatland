%NEE repitition generator for Monte Carlo simulation for EQTP Peatland
%J. Chi
%Apr 21, 2021

clear all; clc; warning off;
tic;
disp 'importing non-gapfilled data...';
mydat = readtable('C:\Users\jihi0001\OneDrive - Sveriges Lantbruksuniversitet\Publications\Peng_Chi_Manuscript\JGR-A_submission\GHG_EQTP_final.csv');
CO2Flux = mydat.co2_flux_filtered_umolm2s;
disp 'done.';

disp 'starting random number generator...';
for i = 1: 100 %update the number of the repitition as needed, here 100
mu = 0;
for j = 1:length(CO2Flux)
    if CO2Flux(j,1) >= 0
        sigma(j,1) = 0.38 + 0.3*CO2Flux(j,1); 
        y(j,1) = laprnd(1, 1, mu, sigma(j,1)); %double exponential distrubtion of CO2 flux residuals from Richardson and Hollinger 2007
        CO2Flux_noise(j,i) = CO2Flux(j,1) + y(j,1);
    else
        sigma(j,1) = 0.47 - 0.12*CO2Flux(j,1);
        y(j,1) = laprnd(1, 1, mu, sigma(j,1));
        CO2Flux_noise(j,i) = CO2Flux(j,1) +y (j,1);
    end
end
end
disp 'done!';

dy=day(mydat.TIMESTAMP,'dayofyear');
disp 'writing output files...';
for i = 1:100 %update the number of the repitition as needed, here 100
    i
    txtfile = ['C:\Users\jihi0001\OneDrive - Sveriges Lantbruksuniversitet\Publications\Peng_Chi_Manuscript\Data\MonteCarlo_input\EQTP_InputData_sim' num2str(i) '.txt']; %update the file path and file name
    f= fopen(txtfile, 'w'); %update output file name
    fprintf(f, 'Year\tDoY\tHour\tNEE\tRg\tTair\tVPD\n  --\t--\t--\tumolm-2s-1\tWm-2\tdegC\thPa\n'); 
    for a = 1: height(mydat)
        fprintf(f, '%f	%f	%f	%f	%f	%f	%f\n',...
            mydat.TIMESTAMP.Year(a,1), dy(a,1), mydat.TIMESTAMP.Hour(a,1)+mydat.TIMESTAMP.Minute(a,1)/60,...
            CO2Flux_noise(a,i), mydat.Rg_Wm2(a,1), mydat.Ta_degC(a,1), mydat.VPD_kPa(a,1)*10);
    end
    fclose(f);
end
disp 'done.';
toc;


%%%%%%%%%%%Uncertainty Analysis (run after gapfilling)%%%%%%%%%%%%%%%%
clear all; clc; close all; warning off;
DataPath1 = 'C:\Users\jihi0001\OneDrive - Sveriges Lantbruksuniversitet\Publications\Peng_Chi_Manuscript\Data\MonteCarlo_output\'; %update the processig file directory
d1 = dir([DataPath1 '*.txt']); 
hd1 = readtable(fullfile(DataPath1, d1(1).name), 'delimiter', 'tab');

t1 = datetime(2013, 12, 1, 0, 0, 0); %update the starting time, begining of the averaging period
t2 = datetime(2016, 7, 31, 23, 30, 0);
ts = datevec((t1: minutes(30): t2)');
[a2, b2, b2] = unique(ts(:, 1: 2), 'rows');

for n = 1: numel(d1)
    dat1 = readtable(fullfile(DataPath1, d1(n).name), 'Delimiter', 'tab', 'HeaderLines', 1);
    dat1.Properties.VariableNames = hd1.Properties.VariableNames;
    dat1 = standardizeMissing(dat1, -9999);
    nee(:, n) = dat1.NEE_f; %100 repetitions of nee
    gpp(:, n) = dat1.GPP_f;
    reco(:, n) = dat1.Reco;
    nee_mon(:, n) = accumarray(b2,dat1.NEE_f,[],@sum)*12*1800/1000000; %100 repetitions of monthly nee sums
    gpp_mon(:, n) = accumarray(b2,dat1.GPP_f,[],@sum)*12*1800/1000000;
    reco_mon(:, n) = accumarray(b2,dat1.Reco,[],@sum)*12*1800/1000000;
end

nee_mon_std = std(nee_mon'); %standard deviation of the 100 monthly nee sums
gpp_mon_std = std(gpp_mon');
reco_mon_std = std(gpp_mon');

nee_std=std(sum(nee)*12*1800/1000000) %standard deviation of the 100 annual nee sums
gpp_std=std(sum(gpp)*12*1800/1000000)
reco_std=std(sum(reco)*12*1800/1000000)
