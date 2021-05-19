%Eastern TQP GHG study
%J. Chi
%Jul 15, 2020
%Updated May 10, 2021

clear all; clc; close all;
warning off;
mydat0 = readtable('C:/Users/jihi0001/Dropbox/Publications/Peng_Chi_Manuscript/Data/GHG_EQTP_final.csv');
mydat = table2timetable(mydat0);
t1 = datetime(2013, 12, 1, 0, 0, 0);
t2 = datetime(2016, 7, 31, 23, 30, 0);
mydat.TIMESTAMP = datetime(t1:minutes(30):t2)';
mydat.rand_err_co2(isnan(mydat.co2_flux_filtered_umolm2s))=0;
mydat.rand_err_ch4(isnan(mydat.ch4_flux_filtered_umolm2s))=0;
mydat.VPD_kPa = mydat.VPD_hPa*0.1;

%Uncertainy analysis due to random measurement errors and gapfilling
ch4_be = mydat.ch4_flux_modeled_umolms - mydat.ch4_flux_filled_umolm2s;
mydat.co2_un = sqrt(mydat.rand_err_co2.^2 + mydat.co2_gapfill_uncertainty.^2);
mydat.ch4_un = sqrt(mydat.rand_err_ch4.^2 + ch4_be.^2);

%30-min non-gapfilled data
mydat.co2 = mydat.co2_flux_filtered_umolm2s*44/1000; %mg CO2-eq m-2 s-1
mydat.ch4_eq1 = mydat.ch4_flux_filtered_umolm2s*16*28/1000; %mg CO2-eq m-2 s-1 for GWP28
mydat.ch4_eq2 = mydat.ch4_flux_filtered_umolm2s*16*45/1000; % mg CO2-eq m-2 s-1 for SGWP45
mydat.ghg1 = mydat.co2+mydat.ch4_eq1;
mydat.ghg2 = mydat.co2+mydat.ch4_eq2;


figure(1)
hold on;
plot(mydat.TIMESTAMP, mydat.co2_flux_filled_umolm2s, '.');
plot(mydat.TIMESTAMP, mydat.co2_flux_filtered_umolm2s, '.');


figure(2)
hold on;
plot(mydat.TIMESTAMP, mydat.ch4_flux_filled_umolm2s, '.');
plot(mydat.TIMESTAMP, mydat.ch4_flux_filtered_umolm2s, '.');


%daily mean or sum data
DailyMean = retime(mydat, 'daily', 'mean');
DailySum = retime(mydat, 'daily', 'sum');

co2_d = DailySum.co2_flux_filled_umolm2s*44*1800/1000000; %unit in g co2 m-2 d-1
ch4_d1 = DailySum.ch4_flux_filled_umolm2s*16*28*1800/1000000; %unit in g co2-eq m-2 d-1 GWP28
ch4_d2 = DailySum.ch4_flux_filled_umolm2s*16*45*1800/1000000; %unit in g co2-eq m-2 d-1 SWGP45
ghg_d1 = co2_d+ch4_d1;
ghg_d2 = co2_d+ch4_d2;


%monthly mean or sum data
MonMean = retime(mydat, 'monthly', 'mean'); %mean
MonSum = retime(mydat, 'monthly', 'sum'); %sum
MonSumsqr = retime(mydat, 'monthly', @sumsqr); %sum square for uncertainty

co2_monthly = MonSum.co2_flux_filled_umolm2s*44*1800/1000000; %unit in g co2-eq m-2 month-1
co2_monthly_sd = sqrt(MonSumsqr.co2_un)*44*1800/1000000; %unit in g co2-eq m-2 month-1
ch4_monthly1 = MonSum.ch4_flux_filled_umolm2s*16*28*1800/1000000; %unit in g co2-eq m-2 month-1
ch4_monthly1_sd = sqrt(MonSumsqr.ch4_un)*16*28*1800/1000000; %unit in g co2-eq m-2 month-1
ch4_monthly2 = MonSum.ch4_flux_filled_umolm2s*16*45*1800/1000000; 
ch4_monthly2_sd =sqrt(MonSumsqr.ch4_un)*16*45*1800/1000000;
ghg_monthly1 = co2_monthly+ch4_monthly1;
ghg_monthly1_sd = sqrt(co2_monthly_sd.^2 + ch4_monthly1_sd.^2);
ghg_monthly2 = co2_monthly+ch4_monthly2;
ghg_monthly2_sd = sqrt(co2_monthly_sd.^2 + ch4_monthly2_sd.^2);


figure(3) %monthly met plots
set(gcf, 'Position', [100 100 700 800]);
subplot(3,1,1);
hold on; box on;
yyaxis left
h1 = bar(MonMean.Rg_Wm2, 'FaceColor', [1.00, 0.41, 0.16], 'BarWidth', 0.6);
xline(1.5, '--k'); xline(13.5, '--k'); xline(25.5, '--k');
ylim([0 400]);
ylabel('R_g (W m^{-2})');
set(gca, 'XTick', 1:32, 'XTickLabel', {'D', 'J', 'F', 'M', 'A', 'M', 'J', 'J',...
    'A', 'S', 'O', 'N', 'D', 'J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D',...
    'J', 'F', 'M', 'A', 'M', 'J', 'J', 'A'},...
    'YTick', [0 100 200 300 400]);
annotation('textbox', [0.27, 0.92, 0.08, 0.04], 'String', "2014",...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle',...
    'FontName', 'Times New Roman', 'FontSize', 14, 'FontWeight', 'Bold');
annotation('textbox', [0.548, 0.92, 0.08, 0.04], 'String', "2015",...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle',...
    'FontName', 'Times New Roman', 'FontSize', 14, 'FontWeight', 'Bold');
annotation('textbox', [0.77, 0.92, 0.08, 0.04], 'String', "2016",...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle',...
    'FontName', 'Times New Roman', 'FontSize', 14, 'FontWeight', 'Bold');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14, 'Layer', 'top');
yyaxis right
h2 = plot(MonMean.Ta_degC, '-o', 'Color', [0.64, 0.08, 0.18], 'LineWidth', 1.5, 'MarkerFaceColor', [0.64, 0.08, 0.18]);
h3 = plot(0:33, zeros(size(0:33)), 'k:');
ylim([-40 20]);
ylabel('T_{air} (\circC)');
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14, 'Layer', 'top');

subplot(3,1,2);
hold on; box on;
yyaxis right
plot(MonMean.VPD_hPa/10, '-o', 'LineWidth', 1.5, 'Color', [0.00, 0.00, 1.00], 'MarkerFaceColor', [0.00, 0.00, 1.00]);
xline(1.5, '--k'); xline(13.5, '--k'); xline(25.5, '--k');
ylim([0 .5]);
ylabel('VPD (kPa)');
set(gca, 'XTick', 1:32, 'XTickLabel', {'D', 'J', 'F', 'M', 'A', 'M', 'J', 'J',...
    'A', 'S', 'O', 'N', 'D', 'J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D',...
    'J', 'F', 'M', 'A', 'M', 'J', 'J', 'A'}, 'YTick', 0:.1:.5);
annotation('textbox', [0.27, 0.322, 0.08, 0.04], 'String', "2014",...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle',...
    'FontName', 'Times New Roman', 'FontSize', 14, 'FontWeight', 'Bold');
annotation('textbox', [0.548, 0.322, 0.08, 0.04], 'String', "2015",...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle',...
    'FontName', 'Times New Roman', 'FontSize', 14, 'FontWeight', 'Bold');
annotation('textbox', [0.77, 0.322, 0.08, 0.04], 'String', "2016",...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle',...
    'FontName', 'Times New Roman', 'FontSize', 14, 'FontWeight', 'Bold');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14, 'Layer', 'top');
yyaxis left
bar(MonSum.PPT_mm, 'FaceColor', [0.30, 0.75, 0.93], 'BarWidth', 0.6);
ylim([0 300]);
ylabel('PPT (mm)');
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
set(gca, 'YTick', 0:50:300);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14, 'Layer', 'top');

subplot(3,1,3);
hold on; box on;
yyaxis right
h1 = plot(MonMean.Soil_10_degC, '-o', 'LineWidth', 1.5, 'Color', [0.35,0.50,0.13], 'MarkerFaceColor', [0.35, 0.50, 0.13]);
h2 = plot(MonMean.Soil_25_degC, '-o', 'LineWidth', 1.5, 'Color', [0.50,0.50,0.50], 'MarkerFaceColor', [0.50,0.50,0.50]);
h3 = plot(MonMean.Soil_40_degC, '-o', 'LineWidth', 1.5, 'Color', 'k', 'MarkerFaceColor', 'k');
h4 = plot(0:33, zeros(size(0:33)), 'k:');
xline(1.5, '--k'); xline(13.5, '--k'); xline(25.5, '--k');
ylim([-5 20]);
set(gca, 'XTick', 1:32, 'XTickLabel', {'D', 'J', 'F', 'M', 'A', 'M', 'J', 'J',...
    'A', 'S', 'O', 'N', 'D', 'J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D',...
    'J', 'F', 'M', 'A', 'M', 'J', 'J', 'A'},...
    'YTick', [-5 0 5 10 15 20]);
annotation('textbox', [0.27, 0.62, 0.08, 0.04], 'String', "2014",...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle',...
    'FontName', 'Times New Roman', 'FontSize', 14, 'FontWeight', 'Bold');
annotation('textbox', [0.548, 0.62, 0.08, 0.04], 'String', "2015",...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle',...
    'FontName', 'Times New Roman', 'FontSize', 14, 'FontWeight', 'Bold');
annotation('textbox', [0.77, 0.62, 0.08, 0.04], 'String', "2016",...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle',...
    'FontName', 'Times New Roman', 'FontSize', 14, 'FontWeight', 'Bold');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14, 'Layer', 'top');
ylabel('T_{soil} (\circC)');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14, 'Layer', 'top');
yyaxis left
bar(MonMean.SWC_10_m3m3, 'FaceColor', [0.64, 0.88, 0.46], 'BarWidth', 0.6);
ylim([0 1]);
ylabel('SWC 10 cm (m^{3} m^{-3})');
legend([h1 h2 h3], {'10 cm', '25 cm', '40 cm'}, 'Location', 'northwest', 'FontSize', 12);
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
set(gca, 'YTick', 0:0.2:1);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14, 'Layer', 'top');


figure(4)
set(gcf, 'Position', [100 100 1000 450]);
hold on; box on;
model_series = [co2_monthly, ch4_monthly1];
model_error = [co2_monthly_sd, ch4_monthly1_sd];
b = bar(model_series, 'grouped', 'EdgeColor', 'none');
b(1).FaceColor = [0.93, 0.69, 0.13];
[ngroups, nbars] = size(model_series);
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, model_series(:,i), model_error(:,i), 'k', 'linestyle', 'none', 'linewidth', 1);
end
h = plot(ghg_monthly1, 'k-o', 'MarkerFaceColor', 'k', 'MarkerSize', 5);
errorbar(ghg_monthly1, ghg_monthly1_sd, 'k', 'linestyle', 'none', 'linewidth', 1);
xline(1.5, '--k');
xline(13.5, '--k');
xline(25.5, '--k');
ylim([-400 300]);
ylabel('CO_2-eq flux (g CO_2-eq m^{-2} month^{-1})');
legend([b, h], 'CO_2', 'CH_4', 'net CO_2-eq', 'Location', 'southeast', 'FontSize', 11);
annotation('textbox', [0.259, 0.93, 0.08, 0.04], 'String', "2014",...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle',...
    'FontName', 'Times New Roman', 'FontSize', 14, 'FontWeight', 'Bold');
annotation('textbox', [0.55, 0.93, 0.08, 0.04], 'String', "2015",...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle',...
    'FontName', 'Times New Roman', 'FontSize', 14, 'FontWeight', 'Bold');
annotation('textbox', [0.78, 0.93, 0.08, 0.04], 'String', "2016",...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle',...
    'FontName', 'Times New Roman', 'FontSize', 14, 'FontWeight', 'Bold');
set(gca, 'YTick', -400:100:300, 'XTick', 1:32, 'XTickLabel', {'D', 'J', 'F', 'M', 'A', 'M', 'J', 'J',...
    'A', 'S', 'O', 'N', 'D', 'J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D',...
    'J', 'F', 'M', 'A', 'M', 'J', 'J', 'A'},...
    'FontName', 'Times New Roman', 'FontSize', 14, 'Layer', 'top');


figure(5)
set(gcf, 'Position', [100 100 1000 450]);
hold on; box on;
model_series = [co2_monthly, ch4_monthly2];
model_error = [co2_monthly_sd, ch4_monthly2_sd];
b = bar(model_series, 'grouped', 'EdgeColor', 'none');
b(1).FaceColor = [0.93, 0.69, 0.13];
[ngroups, nbars] = size(model_series);
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, model_series(:,i), model_error(:,i), 'k', 'linestyle', 'none', 'linewidth', 1);
end
h = plot(ghg_monthly2, 'k-o', 'MarkerFaceColor', 'k', 'MarkerSize', 5);
errorbar(ghg_monthly2, ghg_monthly2_sd, 'k', 'linestyle', 'none', 'linewidth', 1);
xline(1.5, '--k');
xline(13.5, '--k');
xline(25.5, '--k');
ylim([-500 500]);
ylabel('CO_2-eq flux (g CO_2-eq m^{-2} month^{-1})');
legend([b, h], 'CO_2', 'CH_4', 'net CO_2-eq', 'Location', 'southeast', 'FontSize', 11);
annotation('textbox', [0.259, 0.93, 0.08, 0.04], 'String', "2014",...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle',...
    'FontName', 'Times New Roman', 'FontSize', 14, 'FontWeight', 'Bold');
annotation('textbox', [0.55, 0.93, 0.08, 0.04], 'String', "2015",...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle',...
    'FontName', 'Times New Roman', 'FontSize', 14, 'FontWeight', 'Bold');
annotation('textbox', [0.78, 0.93, 0.08, 0.04], 'String', "2016",...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle',...
    'FontName', 'Times New Roman', 'FontSize', 14, 'FontWeight', 'Bold');
set(gca, 'YTick', -500:100:500, 'XTick', 1:32, 'XTickLabel', {'D', 'J', 'F', 'M', 'A', 'M', 'J', 'J',...
    'A', 'S', 'O', 'N', 'D', 'J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D',...
    'J', 'F', 'M', 'A', 'M', 'J', 'J', 'A'},...
    'FontName', 'Times New Roman', 'FontSize', 14, 'Layer', 'top');


%30min non-gapfilled data for PCA
pca1 = table(mydat.ghg1, mydat.co2, mydat.ch4_eq1,...
    mydat.Ta_degC, mydat.Soil_10_degC, mydat.SWC_10_m3m3, mydat.Rg_Wm2, mydat.VPD_kPa, mydat.PPT_mm, ...
    'VariableNames', {'GHG', 'CO2', 'CH4', 'Ta', 'Ts10', 'SWC10', 'Rg', 'VPD', 'PPT'});
pca1 = rmmissing(pca1);
pca_norm1 = zscore(table2array(pca1));
[coeff1, ~, latent1, ~, explained1] = pca(pca_norm1);
[R1, p1] = corrcoef(pca_norm1);

%daily gapfilled data for PCA
pca2 = table(ghg_d1, co2_d, ch4_d1, DailyMean.Ta_degC, DailyMean.Soil_10_degC, DailyMean.SWC_10_m3m3,...
    DailyMean.Rg_Wm2, DailyMean.VPD_kPa, DailySum.PPT_mm,...
    'VariableNames', {'GHG', 'CO2', 'CH4', 'Ta', 'Ts10', 'SWC10', 'Rg', 'VPD', 'PPT'});
pca2 = rmmissing(pca2);
pca_norm2 = zscore(table2array(pca2));
[coeff2, ~, latent2, ~, explained2] = pca(pca_norm2);
[R2, p2] = corrcoef(pca_norm2);

%monthly gapfilled data for PCA
pca3 = table(ghg_monthly1, co2_monthly, ch4_monthly1,...
    MonMean.Ta_degC, MonMean.Soil_10_degC, MonMean.SWC_10_m3m3,...
    MonMean.Rg_Wm2, MonMean.VPD_kPa, MonSum.PPT_mm,...
    'VariableNames', {'GHG', 'CO2', 'CH4', 'Ta', 'Ts10', 'SWC10', 'Rg', 'VPD', 'PPT'});
pca3 = rmmissing(pca3);
pca_norm3 = zscore(table2array(pca3));
[coeff3, ~, latent3, ~, explained3] = pca(pca_norm3);
[R3, p3] = corrcoef(pca_norm3);

pca4 = table(mydat.ghg2, mydat.co2, mydat.ch4_eq2,...
    mydat.Ta_degC, mydat.Soil_10_degC, mydat.SWC_10_m3m3, mydat.Rg_Wm2, mydat.VPD_kPa, mydat.PPT_mm, ...
    'VariableNames', {'GHG', 'CO2', 'CH4', 'Ta', 'Ts10', 'SWC10', 'Rg', 'VPD', 'PPT'});
pca4 = rmmissing(pca4);
pca_norm4 = zscore(table2array(pca4));
[coeff4, ~, latent4, ~, explained4] = pca(pca_norm4);
[R4, p4] = corrcoef(pca_norm4);

%daily gapfilled data for PCA
pca5 = table(ghg_d2, co2_d, ch4_d2, DailyMean.Ta_degC, DailyMean.Soil_10_degC, DailyMean.SWC_10_m3m3,...
    DailyMean.Rg_Wm2, DailyMean.VPD_kPa, DailySum.PPT_mm,...
    'VariableNames', {'GHG', 'CO2', 'CH4', 'Ta', 'Ts10', 'SWC10', 'Rg', 'VPD', 'PPT'});
pca5 = rmmissing(pca5);
pca_norm5 = zscore(table2array(pca5));
[coeff5, ~, latent5, ~, explained5] = pca(pca_norm5);
[R5, p5] = corrcoef(pca_norm5);

%monthly gapfilled data for PCA
pca6 = table(ghg_monthly2, co2_monthly, ch4_monthly2,...
    MonMean.Ta_degC, MonMean.Soil_10_degC, MonMean.SWC_10_m3m3,...
    MonMean.Rg_Wm2, MonMean.VPD_kPa, MonSum.PPT_mm,...
    'VariableNames', {'GHG', 'CO2', 'CH4', 'Ta', 'Ts10', 'SWC10', 'Rg', 'VPD', 'PPT'});
pca6 = rmmissing(pca6);
pca_norm6 = zscore(table2array(pca6));
[coeff6, ~, latent6, ~, explained6] = pca(pca_norm6);
[R6, p6] = corrcoef(pca_norm6);


[X,Y] = meshgrid(1:9, 1:9);
mycolormap = customcolormap(linspace(0,1,2), [0.93 0.98 1; 0 0.45 0.74]); %blue gradient

figure(6)
set(gcf, 'Position', [100 100 1400 820]);
subplot(2,3,1);
h1 = biplot(coeff1(:, 1:2), 'LineWidth', 1.5, 'VarLabels', {'GHG', 'CO_2', 'CH_4', 'T_{air}', 'T_{soil}', 'SWC', 'R_g', 'VPD', 'PPT'});
h1(end).XData = [-0.8 0.8 NaN 0 0];
h1(end).YData = [0 0 NaN -0.8 0.8];
co=[0.00, 0.00, 0.00; 0.93, 0.69, 0.13; 0.85, 0.33, 0.10;...
    0.00, 0.45, 0.74; 0.47, 0.67, 0.19; 0.30, 0.75, 0.93;...
    0.64, 0.08, 0.18; 0.72, 0.27, 1.00; 0.00,0.00,1.00];
for i = 1:9
    set(h1(i), 'Color', co(i, :));
    set(h1(i+9), 'Color', co(i, :));
    set(h1(i+18), 'FontName', 'Times New Roman', 'Color', co(i, :), 'Interpreter', 'tex');
end
xlim([-0.8 0.8]);
ylim([-0.8 0.8]);
xlabel(['PC1 (' num2str(explained1(1),'%02.f') '%)']);
ylabel(['PC2 (' num2str(explained1(2),'%02.f') '%)']);
title('(a) Half-hourly', 'FontSize', 14);
set(gca, 'XTick', -0.8:0.2:0.8, 'YTick', -0.8:0.2:0.8);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
box on;

subplot(2,3,2);
h2 = biplot(coeff2(:, 1:2), 'LineWidth', 1.5, 'VarLabels', {'GHG', 'CO_2', 'CH_4', 'T_{air}', 'T_{soil}', 'SWC', 'R_g', 'VPD', 'PPT'});
h2(end).XData = [-0.8 0.8 NaN 0 0];
h2(end).YData = [0 0 NaN -0.8 0.8];
for i = 1:9
    set(h2(i), 'Color', co(i ,:));
    set(h2(i+9), 'Color', co(i, :));
    set(h2(i+18), 'FontName', 'Times New Roman', 'Color', co(i, :), 'Interpreter', 'tex');
end
xlim([-0.8 0.8]);
ylim([-0.8 0.8]);
xlabel(['PC1 (' num2str(explained2(1),'%02.f') '%)']);
ylabel(['PC2 (' num2str(explained2(2),'%02.f') '%)']);
title('(b) Daily', 'FontSize', 14);
set(gca, 'XTick', -0.8:0.2:0.8, 'YTick', -0.8:0.2:0.8);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
box on;

subplot(2,3,3);
h3 = biplot(coeff3(:, 1:2), 'LineWidth', 1.5, 'VarLabels', {'GHG', 'CO_2', 'CH_4', 'T_{air}', 'T_{soil}', 'SWC', 'R_g', 'VPD', 'PPT'});
h3(end).XData = [-0.8 0.8 NaN 0 0];
h3(end).YData = [0 0 NaN -0.8 0.8];
for i = 1:9
    set(h3(i), 'Color', co(i, :));
    set(h3(i+9), 'Color', co(i, :));
    set(h3(i+18), 'FontName', 'Times New Roman', 'Color', co(i, :), 'Interpreter', 'tex');
end
xlim([-0.8 0.8]);
ylim([-0.8 0.8]);
xlabel(['PC1 (' num2str(explained3(1),'%02.f') '%)']);
ylabel(['PC2 (' num2str(explained3(2),'%02.f') '%)']);
title('(c) Monthly', 'FontSize', 14);
set(gca, 'XTick', -0.8:0.2:0.8, 'YTick', -0.8:0.2:0.8);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
box on;

subplot(2,3,4);
h1 = imagesc(X(:), Y(:), R1);
colormap(mycolormap);
caxis([-1 1]);
txt1 = sprintfc('%.2f', R1);
for i = 1:9
    for j = 1:9
        if p1(i,j) < 0.001
            text(X(i,j), Y(i,j), txt1{i,j}, 'horizontalalignment', 'center', 'verticalalignment', 'middle',...
                'Color', 'r', 'FontName', 'Times New Roman');
        else
            text(X(i,j), Y(i,j), txt1{i,j}, 'horizontalalignment', 'center', 'verticalalignment', 'middle',...
                'Color', 'k', 'FontName', 'Times New Roman');
        end
    end
end
title('(d) Half-hourly', 'FontSize', 14);
set(gca, 'XTick', 1:9, 'XTickLabel', {'GHG', 'CO_2', 'CH_4', 'T_{air}', 'T_{soil}', 'SWC', 'R_g', 'VPD', 'PPT'},...
    'YTick', 1:9, 'YTickLabel', {'GHG', 'CO_2', 'CH_4', 'T_{air}', 'T_{soil}', 'SWC', 'R_g', 'VPD', 'PPT'},...
    'FontName', 'Times New Roman');

subplot(2,3,5);
h2 = imagesc(X(:), Y(:), R2);
colormap(mycolormap);
caxis([-1 1]);
txt2 = sprintfc('%.2f', R2);
for i = 1:9
    for j = 1:9
        if p2(i,j) < 0.001
            text(X(i,j), Y(i,j), txt2{i,j}, 'horizontalalignment', 'center', 'verticalalignment', 'middle',...
                'Color', 'r', 'FontName', 'Times New Roman');
        else
            text(X(i,j), Y(i,j), txt2{i,j}, 'horizontalalignment', 'center', 'verticalalignment', 'middle',...
                'Color', 'k', 'FontName', 'Times New Roman');
        end
    end
end
title('(e) Daily', 'FontSize', 14);
set(gca, 'XTick', 1:9, 'XTickLabel', {'GHG', 'CO_2', 'CH_4', 'T_{air}', 'T_{soil}', 'SWC', 'R_g', 'VPD', 'PPT'},...
    'YTick', 1:9, 'YTickLabel', {'GHG', 'CO_2', 'CH_4', 'T_{air}', 'T_{soil}', 'SWC', 'R_g', 'VPD', 'PPT'},...
    'FontName', 'Times New Roman');

subplot(2,3,6);
h3 = imagesc(X(:), Y(:), R3);
colormap(mycolormap);
caxis([-1 1]);
colorbar
txt3 = sprintfc('%.2f', R3);
for i = 1:9
    for j = 1:9
        if p3(i,j) < 0.001
            text(X(i,j), Y(i,j), txt3{i,j}, 'horizontalalignment', 'center', 'verticalalignment', 'middle',...
                'Color', 'r', 'FontName', 'Times New Roman');
        else
            text(X(i,j), Y(i,j), txt3{i,j}, 'horizontalalignment', 'center', 'verticalalignment', 'middle',...
                'Color', 'k', 'FontName', 'Times New Roman');
        end
    end
end
title('(f) Monthly', 'FontSize', 14);
set(gca, 'XTick', 1:9, 'XTickLabel', {'GHG', 'CO_2', 'CH_4', 'T_{air}', 'T_{soil}', 'SWC', 'R_g', 'VPD', 'PPT'},...
    'YTick', 1:9, 'YTickLabel', {'GHG', 'CO_2', 'CH_4', 'T_{air}', 'T_{soil}', 'SWC', 'R_g', 'VPD', 'PPT'},...
    'FontName', 'Times New Roman');


figure(7)
set(gcf, 'Position', [100 100 1400 820]);
subplot(2,3,1);
h1 = biplot(coeff4(:, 1:2), 'LineWidth', 1.5, 'VarLabels', {'GHG', 'CO_2', 'CH_4', 'T_{air}', 'T_{soil}', 'SWC', 'R_g', 'VPD', 'PPT'});
h1(end).XData = [-0.8 0.8 NaN 0 0];
h1(end).YData = [0 0 NaN -0.8 0.8];
co=[0.00, 0.00, 0.00; 0.93, 0.69, 0.13; 0.85, 0.33, 0.10;...
    0.00, 0.45, 0.74; 0.47, 0.67, 0.19; 0.30, 0.75, 0.93;...
    0.64, 0.08, 0.18; 0.72, 0.27, 1.00; 0.00,0.00,1.00];
for i = 1:9
    set(h1(i), 'Color', co(i, :));
    set(h1(i+9), 'Color', co(i, :));
    set(h1(i+18), 'FontName', 'Times New Roman', 'Color', co(i, :), 'Interpreter', 'tex');
end
xlim([-0.8 0.8]);
ylim([-0.8 0.8]);
xlabel(['PC1 (' num2str(explained4(1),'%02.f') '%)']);
ylabel(['PC2 (' num2str(explained4(2),'%02.f') '%)']);
title('(a) Half-hourly', 'FontSize', 14);
set(gca, 'XTick', -0.8:0.2:0.8, 'YTick', -0.8:0.2:0.8);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
box on;

subplot(2,3,2);
h2 = biplot(coeff5(:, 1:2), 'LineWidth', 1.5, 'VarLabels', {'GHG', 'CO_2', 'CH_4', 'T_{air}', 'T_{soil}', 'SWC', 'R_g', 'VPD', 'PPT'});
h2(end).XData = [-0.8 0.8 NaN 0 0];
h2(end).YData = [0 0 NaN -0.8 0.8];
for i = 1:9
    set(h2(i), 'Color', co(i ,:));
    set(h2(i+9), 'Color', co(i, :));
    set(h2(i+18), 'FontName', 'Times New Roman', 'Color', co(i, :), 'Interpreter', 'tex');
end
xlim([-0.8 0.8]);
ylim([-0.8 0.8]);
xlabel(['PC1 (' num2str(explained5(1),'%02.f') '%)']);
ylabel(['PC2 (' num2str(explained5(2),'%02.f') '%)']);
title('(b) Daily', 'FontSize', 14);
set(gca, 'XTick', -0.8:0.2:0.8, 'YTick', -0.8:0.2:0.8);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
box on;

subplot(2,3,3);
h3 = biplot(coeff6(:, 1:2), 'LineWidth', 1.5, 'VarLabels', {'GHG', 'CO_2', 'CH_4', 'T_{air}', 'T_{soil}', 'SWC', 'R_g', 'VPD', 'PPT'});
h3(end).XData = [-0.8 0.8 NaN 0 0];
h3(end).YData = [0 0 NaN -0.8 0.8];
for i = 1:9
    set(h3(i), 'Color', co(i, :));
    set(h3(i+9), 'Color', co(i, :));
    set(h3(i+18), 'FontName', 'Times New Roman', 'Color', co(i, :), 'Interpreter', 'tex');
end
xlim([-0.8 0.8]);
ylim([-0.8 0.8]);
xlabel(['PC1 (' num2str(explained6(1),'%02.f') '%)']);
ylabel(['PC2 (' num2str(explained6(2),'%02.f') '%)']);
title('(c) Monthly', 'FontSize', 14);
set(gca, 'XTick', -0.8:0.2:0.8, 'YTick', -0.8:0.2:0.8);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
box on;

subplot(2,3,4);
h1 = imagesc(X(:), Y(:), R4);
colormap(mycolormap);
caxis([-1 1]);
txt4 = sprintfc('%.2f', R4);
for i = 1:9
    for j = 1:9
        if p4(i,j) < 0.001
            text(X(i,j), Y(i,j), txt4{i,j}, 'horizontalalignment', 'center', 'verticalalignment', 'middle',...
                'Color', 'r', 'FontName', 'Times New Roman');
        else
            text(X(i,j), Y(i,j), txt4{i,j}, 'horizontalalignment', 'center', 'verticalalignment', 'middle',...
                'Color', 'k', 'FontName', 'Times New Roman');
        end
    end
end
title('(d) Half-hourly', 'FontSize', 14);
set(gca, 'XTick', 1:9, 'XTickLabel', {'GHG', 'CO_2', 'CH_4', 'T_{air}', 'T_{soil}', 'SWC', 'R_g', 'VPD', 'PPT'},...
    'YTick', 1:9, 'YTickLabel', {'GHG', 'CO_2', 'CH_4', 'T_{air}', 'T_{soil}', 'SWC', 'R_g', 'VPD', 'PPT'},...
    'FontName', 'Times New Roman');

subplot(2,3,5);
h2 = imagesc(X(:), Y(:), R5);
colormap(mycolormap);
caxis([-1 1]);
txt5 = sprintfc('%.2f', R5);
for i = 1:9
    for j = 1:9
        if p5(i,j) < 0.001
            text(X(i,j), Y(i,j), txt5{i,j}, 'horizontalalignment', 'center', 'verticalalignment', 'middle',...
                'Color', 'r', 'FontName', 'Times New Roman');
        else
            text(X(i,j), Y(i,j), txt5{i,j}, 'horizontalalignment', 'center', 'verticalalignment', 'middle',...
                'Color', 'k', 'FontName', 'Times New Roman');
        end
    end
end
title('(e) Daily', 'FontSize', 14);
set(gca, 'XTick', 1:9, 'XTickLabel', {'GHG', 'CO_2', 'CH_4', 'T_{air}', 'T_{soil}', 'SWC', 'R_g', 'VPD', 'PPT'},...
    'YTick', 1:9, 'YTickLabel', {'GHG', 'CO_2', 'CH_4', 'T_{air}', 'T_{soil}', 'SWC', 'R_g', 'VPD', 'PPT'},...
    'FontName', 'Times New Roman');

subplot(2,3,6);
h3 = imagesc(X(:), Y(:), R6);
colormap(mycolormap);
caxis([-1 1]);
colorbar
txt6 = sprintfc('%.2f', R6);
for i = 1:9
    for j = 1:9
        if p6(i,j) < 0.001
            text(X(i,j), Y(i,j), txt6{i,j}, 'horizontalalignment', 'center', 'verticalalignment', 'middle',...
                'Color', 'r', 'FontName', 'Times New Roman');
        else
            text(X(i,j), Y(i,j), txt3{i,j}, 'horizontalalignment', 'center', 'verticalalignment', 'middle',...
                'Color', 'k', 'FontName', 'Times New Roman');
        end
    end
end
title('(f) Monthly', 'FontSize', 14);
set(gca, 'XTick', 1:9, 'XTickLabel', {'GHG', 'CO_2', 'CH_4', 'T_{air}', 'T_{soil}', 'SWC', 'R_g', 'VPD', 'PPT'},...
    'YTick', 1:9, 'YTickLabel', {'GHG', 'CO_2', 'CH_4', 'T_{air}', 'T_{soil}', 'SWC', 'R_g', 'VPD', 'PPT'},...
    'FontName', 'Times New Roman');

st1_ts = timerange('03/06/2014 00:00', '04/13/2014 00:00');
st1 = mydat(st1_ts, :); %soil thawing 2014

gs1_ts = timerange('04/13/2014 00:00', '10/02/2014 00:00');
gs1 = mydat(gs1_ts, :); %growing season 2014 

sf1_ts = timerange('10/02/2014 00:00', '12/03/2014 00:00');
sf1 = mydat(sf1_ts, :); %soil freezing 2014

w1 = mydat(mydat.TIMESTAMP.Year==2014, :);
w1(st1_ts, :)=[]; w1(gs1_ts, :)=[]; w1(sf1_ts, :)=[];

st2_ts = timerange('03/14/2015 00:00', '05/12/2015 00:00');
st2 = mydat(st2_ts, :);

gs2_ts = timerange('05/12/2015 00:00', '10/08/2015 00:00');
gs2 = mydat(gs2_ts, :);

sf2_ts = timerange('10/08/2015 00:00', '12/17/2015 00:00');
sf2 = mydat(sf2_ts, :);

w2 = mydat(mydat.TIMESTAMP.Year==2015, :);
w2(st2_ts, :)=[]; w2(gs2_ts, :)=[]; w2(sf2_ts, :)=[];

st3_ts = timerange('03/02/2016 00:00', '05/01/2016 00:00');
st3 = mydat(st3_ts, :);

gs3_ts = timerange('05/01/2016 00:00', '08/01/2016 00:00');
gs3 = mydat(gs3_ts, :);

w3 = mydat(mydat.TIMESTAMP.Year==2016, :);
w3(st3_ts, :)=[]; w2(gs3_ts, :)=[];

w0 = mydat(mydat.TIMESTAMP.Year==2013, :);

w = [w0; w1; w2];
st = [st1; st2; st3];
gs = [gs1; gs2; gs3];
sf = [sf1; sf2];

%%%%%%%%%%%%%Rg responses%%%%%%%%%%%%%%%
bn = linspace(0, round(max(w.Rg_Wm2), -2), 16);
for i = 1:15
    bin_Rg1{i} = w.Rg_Wm2(w.Rg_Wm2 >= bn(i) & w.Rg_Wm2 < bn(i+1));
    bin_Rg_avg1(i,1) = nanmean(bin_Rg1{i});
    Rg_co2_1{i} = w.co2(w.Rg_Wm2 >= bn(i) & w.Rg_Wm2 < bn(i+1));
    Rg_ch4_1a{i} = w.ch4_eq1(w.Rg_Wm2 >= bn(i) & w.Rg_Wm2 < bn(i+1));
    Rg_ch4_1b{i} = w.ch4_eq2(w.Rg_Wm2 >= bn(i) & w.Rg_Wm2 < bn(i+1));
    Rg_ghg1a{i} = w.ghg1(w.Rg_Wm2 >= bn(i) & w.Rg_Wm2 < bn(i+1));
    Rg_ghg1b{i} = w.ghg2(w.Rg_Wm2 >= bn(i) & w.Rg_Wm2 < bn(i+1));
    Rg_co2_avg1(i,1) = nanmean(Rg_co2_1{i});
    Rg_ch4_avg1a(i,1) = nanmean(Rg_ch4_1a{i});
    Rg_ch4_avg1b(i,1) = nanmean(Rg_ch4_1b{i});
    Rg_ghg_avg1a(i,1) = nanmean(Rg_ghg1a{i});
    Rg_ghg_avg1b(i,1) = nanmean(Rg_ghg1b{i});
    Rg_co2_std1(i,1) = nanstd(Rg_co2_1{i});
    Rg_ch4_std1a(i,1) = nanstd(Rg_ch4_1a{i});
    Rg_ch4_std1b(i,1) = nanstd(Rg_ch4_1b{i});
    Rg_ghg_std1a(i,1) = nanstd(Rg_ghg1a{i});
    Rg_ghg_std1b(i,1) = nanstd(Rg_ghg1b{i});
    Rg_co2_se1(i,1) = Rg_co2_std1(i,1)/sqrt(length(Rg_co2_1{i}));
    Rg_ch4_se1a(i,1) = Rg_ch4_std1a(i,1)/sqrt(length(Rg_ch4_1a{i}));
    Rg_ch4_se1b(i,1) = Rg_ch4_std1b(i,1)/sqrt(length(Rg_ch4_1b{i}));
    Rg_ghg_se1a(i,1) = Rg_ghg_std1a(i,1)/sqrt(length(Rg_ghg1a{i}));
    Rg_ghg_se1b(i,1) = Rg_ghg_std1b(i,1)/sqrt(length(Rg_ghg1b{i}));
end

bn = linspace(0, round(max(st.Rg_Wm2), -2), 16);
for i = 1:15
    bin_Rg2{i} = st.Rg_Wm2(st.Rg_Wm2 >= bn(i) & st.Rg_Wm2 < bn(i+1));
    bin_Rg_avg2(i,1) = nanmean(bin_Rg2{i});
    Rg_co2_2{i} = st.co2(st.Rg_Wm2 >= bn(i) & st.Rg_Wm2 < bn(i+1));
    Rg_ch4_2a{i} = st.ch4_eq1(st.Rg_Wm2 >= bn(i) & st.Rg_Wm2 < bn(i+1));
    Rg_ch4_2b{i} = st.ch4_eq2(st.Rg_Wm2 >= bn(i) & st.Rg_Wm2 < bn(i+1));
    Rg_ghg2a{i} = st.ghg1(st.Rg_Wm2 >= bn(i) & st.Rg_Wm2 < bn(i+1));
    Rg_ghg2b{i} = st.ghg2(st.Rg_Wm2 >= bn(i) & st.Rg_Wm2 < bn(i+1));
    Rg_co2_avg2(i,1) = nanmean(Rg_co2_2{i});
    Rg_ch4_avg2a(i,1) = nanmean(Rg_ch4_2a{i});
    Rg_ch4_avg2b(i,1) = nanmean(Rg_ch4_2b{i});
    Rg_ghg_avg2a(i,1) = nanmean(Rg_ghg2a{i});
    Rg_ghg_avg2b(i,1) = nanmean(Rg_ghg2b{i});
    Rg_co2_std2(i,1) = nanstd(Rg_co2_2{i});
    Rg_ch4_std2a(i,1) = nanstd(Rg_ch4_2a{i});
    Rg_ch4_std2b(i,1) = nanstd(Rg_ch4_2b{i});
    Rg_ghg_std2a(i,1) = nanstd(Rg_ghg2a{i});
    Rg_ghg_std2b(i,1) = nanstd(Rg_ghg2b{i});
    Rg_co2_se2(i,1) = Rg_co2_std2(i,1)/sqrt(length(Rg_co2_2{i}));
    Rg_ch4_se2a(i,1) = Rg_ch4_std2a(i,1)/sqrt(length(Rg_ch4_2a{i}));
    Rg_ch4_se2b(i,1) = Rg_ch4_std2b(i,1)/sqrt(length(Rg_ch4_2b{i}));
    Rg_ghg_se2a(i,1) = Rg_ghg_std2a(i,1)/sqrt(length(Rg_ghg2a{i}));
    Rg_ghg_se2b(i,1) = Rg_ghg_std2b(i,1)/sqrt(length(Rg_ghg2b{i}));
end

bn = linspace(0, round(max(gs.Rg_Wm2), -2), 16);
for i = 1:15
    bin_Rg3{i} = gs.Rg_Wm2(gs.Rg_Wm2 >= bn(i) & gs.Rg_Wm2 < bn(i+1));
    bin_Rg_avg3(i,1) = nanmean(bin_Rg3{i});
    Rg_co2_3{i} = gs.co2(gs.Rg_Wm2 >= bn(i) & gs.Rg_Wm2 < bn(i+1));
    Rg_ch4_3a{i} = gs.ch4_eq1(gs.Rg_Wm2 >= bn(i) & gs.Rg_Wm2 < bn(i+1));
    Rg_ch4_3b{i} = gs.ch4_eq2(gs.Rg_Wm2 >= bn(i) & gs.Rg_Wm2 < bn(i+1));
    Rg_ghg3a{i} = gs.ghg1(gs.Rg_Wm2 >= bn(i) & gs.Rg_Wm2 < bn(i+1));
    Rg_ghg3b{i} = gs.ghg2(gs.Rg_Wm2 >= bn(i) & gs.Rg_Wm2 < bn(i+1));
    Rg_co2_avg3(i,1) = nanmean(Rg_co2_3{i});
    Rg_ch4_avg3a(i,1) = nanmean(Rg_ch4_3a{i});
    Rg_ch4_avg3b(i,1) = nanmean(Rg_ch4_3b{i});
    Rg_ghg_avg3a(i,1) = nanmean(Rg_ghg3a{i});
    Rg_ghg_avg3b(i,1) = nanmean(Rg_ghg3b{i});
    Rg_co2_std3(i,1) = nanstd(Rg_co2_3{i});
    Rg_ch4_std3a(i,1) = nanstd(Rg_ch4_3a{i});
    Rg_ch4_std3b(i,1) = nanstd(Rg_ch4_3b{i});
    Rg_ghg_std3a(i,1) = nanstd(Rg_ghg3a{i});
    Rg_ghg_std3b(i,1) = nanstd(Rg_ghg3b{i});
    Rg_co2_se3(i,1) = Rg_co2_std3(i,1)/sqrt(length(Rg_co2_3{i}));
    Rg_ch4_se3a(i,1) = Rg_ch4_std3a(i,1)/sqrt(length(Rg_ch4_3a{i}));
    Rg_ch4_se3b(i,1) = Rg_ch4_std3b(i,1)/sqrt(length(Rg_ch4_3b{i}));
    Rg_ghg_se3a(i,1) = Rg_ghg_std3a(i,1)/sqrt(length(Rg_ghg3a{i}));
    Rg_ghg_se3b(i,1) = Rg_ghg_std3b(i,1)/sqrt(length(Rg_ghg3b{i}));
end

bn = linspace(0, round(max(sf.Rg_Wm2), -2), 16);
for i = 1:15
    bin_Rg4{i} = sf.Rg_Wm2(sf.Rg_Wm2 >= bn(i) & sf.Rg_Wm2 < bn(i+1));
    bin_Rg_avg4(i,1) = nanmean(bin_Rg4{i});
    Rg_co2_4{i} = sf.co2(sf.Rg_Wm2 >= bn(i) & sf.Rg_Wm2 < bn(i+1));
    Rg_ch4_4a{i} = sf.ch4_eq1(sf.Rg_Wm2 >= bn(i) & sf.Rg_Wm2 < bn(i+1));
    Rg_ch4_4b{i} = sf.ch4_eq2(sf.Rg_Wm2 >= bn(i) & sf.Rg_Wm2 < bn(i+1));
    Rg_ghg4a{i} = sf.ghg1(sf.Rg_Wm2 >= bn(i) & sf.Rg_Wm2 < bn(i+1));
    Rg_ghg4b{i} = sf.ghg2(sf.Rg_Wm2 >= bn(i) & sf.Rg_Wm2 < bn(i+1));
    Rg_co2_avg4(i,1) = nanmean(Rg_co2_4{i});
    Rg_ch4_avg4a(i,1) = nanmean(Rg_ch4_4a{i});
    Rg_ch4_avg4b(i,1) = nanmean(Rg_ch4_4b{i});
    Rg_ghg_avg4a(i,1) = nanmean(Rg_ghg4a{i});
    Rg_ghg_avg4b(i,1) = nanmean(Rg_ghg4b{i});
    Rg_co2_std4(i,1) = nanstd(Rg_co2_4{i});
    Rg_ch4_std4a(i,1) = nanstd(Rg_ch4_4a{i});
    Rg_ch4_std4b(i,1) = nanstd(Rg_ch4_4b{i});
    Rg_ghg_std4a(i,1) = nanstd(Rg_ghg4a{i});
    Rg_ghg_std4b(i,1) = nanstd(Rg_ghg4b{i});
    Rg_co2_se4(i,1) = Rg_co2_std4(i,1)/sqrt(length(Rg_co2_4{i}));
    Rg_ch4_se4a(i,1) = Rg_ch4_std4a(i,1)/sqrt(length(Rg_ch4_4a{i}));
    Rg_ch4_se4b(i,1) = Rg_ch4_std4b(i,1)/sqrt(length(Rg_ch4_4b{i}));
    Rg_ghg_se4a(i,1) = Rg_ghg_std4a(i,1)/sqrt(length(Rg_ghg4a{i}));
    Rg_ghg_se4b(i,1) = Rg_ghg_std4b(i,1)/sqrt(length(Rg_ghg4b{i}));
end


%%%%%%%%%%%%%%Tsoil response%%%%%%%%%%%%%%%%%%%%
bn = linspace(-5.5, 0.5, 16);
for i = 1:15
    bin_Ts1{i} = w.Soil_10_degC(w.Soil_10_degC >= bn(i) & w.Soil_10_degC < bn(i+1));
    bin_Ts_avg1(i,1) = nanmean(bin_Ts1{i});
    Ts_ghg1a{i} = w.ghg1(w.Soil_10_degC >= bn(i) & w.Soil_10_degC < bn(i+1));
    Ts_ghg1b{i} = w.ghg2(w.Soil_10_degC >= bn(i) & w.Soil_10_degC < bn(i+1));
    Ts_co2_1{i} = w.co2(w.Soil_10_degC >= bn(i) & w.Soil_10_degC < bn(i+1));
    Ts_ch4_1a{i} = w.ch4_eq1(w.Soil_10_degC >= bn(i) & w.Soil_10_degC < bn(i+1));
    Ts_ch4_1b{i} = w.ch4_eq2(w.Soil_10_degC >= bn(i) & w.Soil_10_degC < bn(i+1));
    Ts_ghg_avg1a(i,1) = nanmean(Ts_ghg1a{i});
    Ts_ghg_avg1b(i,1) = nanmean(Ts_ghg1b{i});
    Ts_co2_avg1(i,1) = nanmean(Ts_co2_1{i});
    Ts_ch4_avg1a(i,1) = nanmean(Ts_ch4_1a{i});
    Ts_ch4_avg1b(i,1) = nanmean(Ts_ch4_1b{i});
    Ts_ghg_std1a(i,1) = nanstd(Ts_ghg1a{i});
    Ts_ghg_std1b(i,1) = nanstd(Ts_ghg1b{i});
    Ts_co2_std1(i,1) = nanstd(Ts_co2_1{i});
    Ts_ch4_std1a(i,1) = nanstd(Ts_ch4_1a{i});
    Ts_ch4_std1b(i,1) = nanstd(Ts_ch4_1b{i});
    Ts_ghg_se1a(i,1) = Ts_ghg_std1a(i,1)/sqrt(length(Ts_ghg1a{i}));
    Ts_ghg_se1b(i,1) = Ts_ghg_std1b(i,1)/sqrt(length(Ts_ghg1b{i}));
    Ts_co2_se1(i,1) = Ts_co2_std1(i,1)/sqrt(length(Ts_co2_1{i}));
    Ts_ch4_se1a(i,1) = Ts_ch4_std1a(i,1)/sqrt(length(Ts_ch4_1a{i}));
    Ts_ch4_se1b(i,1) = Ts_ch4_std1b(i,1)/sqrt(length(Ts_ch4_1b{i}));
end

bn = linspace(-1, 6.5, 16);
for i = 1:15
    bin_Ts2{i} = st.Soil_10_degC(st.Soil_10_degC >= bn(i) & st.Soil_10_degC < bn(i+1));
    bin_Ts_avg2(i,1) = nanmean(bin_Ts2{i});
    Ts_ghg2a{i} = st.ghg1(st.Soil_10_degC >= bn(i) & st.Soil_10_degC < bn(i+1));
    Ts_ghg2b{i} = st.ghg2(st.Soil_10_degC >= bn(i) & st.Soil_10_degC < bn(i+1));
    Ts_co2_2{i} = st.co2(st.Soil_10_degC >= bn(i) & st.Soil_10_degC < bn(i+1));
    Ts_ch4_2a{i} = st.ch4_eq1(st.Soil_10_degC >= bn(i) & st.Soil_10_degC < bn(i+1));
    Ts_ch4_2b{i} = st.ch4_eq2(st.Soil_10_degC >= bn(i) & st.Soil_10_degC < bn(i+1));
    Ts_ghg_avg2a(i,1) = nanmean(Ts_ghg2a{i});
    Ts_ghg_avg2b(i,1) = nanmean(Ts_ghg2b{i});
    Ts_co2_avg2(i,1) = nanmean(Ts_co2_2{i});
    Ts_ch4_avg2a(i,1) = nanmean(Ts_ch4_2a{i});
    Ts_ch4_avg2b(i,1) = nanmean(Ts_ch4_2b{i});
    Ts_ghg_std2a(i,1) = nanstd(Ts_ghg2a{i});
    Ts_ghg_std2b(i,1) = nanstd(Ts_ghg2b{i});
    Ts_co2_std2(i,1) = nanstd(Ts_co2_2{i});
    Ts_ch4_std2a(i,1) = nanstd(Ts_ch4_2a{i});
    Ts_ch4_std2b(i,1) = nanstd(Ts_ch4_2b{i});
    Ts_ghg_se2a(i,1) = Ts_ghg_std2a(i,1)/sqrt(length(Ts_ghg2a{i}));
    Ts_ghg_se2b(i,1) = Ts_ghg_std2b(i,1)/sqrt(length(Ts_ghg2b{i}));
    Ts_co2_se2(i,1) = Ts_co2_std2(i,1)/sqrt(length(Ts_co2_2{i}));
    Ts_ch4_se2a(i,1)= Ts_ch4_std2a(i,1)/sqrt(length(Ts_ch4_2a{i}));
    Ts_ch4_se2b(i,1)= Ts_ch4_std2b(i,1)/sqrt(length(Ts_ch4_2b{i}));
end

bn = linspace(0, 23, 16);
for i = 1:15
    bin_Ts3{i} = gs.Soil_10_degC(gs.Soil_10_degC >= bn(i) & gs.Soil_10_degC < bn(i+1));
    bin_Ts_avg3(i,1) = nanmean(bin_Ts3{i});
    Ts_ghg3a{i} = gs.ghg1(gs.Soil_10_degC >= bn(i) & gs.Soil_10_degC < bn(i+1));
    Ts_ghg3b{i} = gs.ghg2(gs.Soil_10_degC >= bn(i) & gs.Soil_10_degC < bn(i+1));
    Ts_co2_3{i} = gs.co2(gs.Soil_10_degC >= bn(i) & gs.Soil_10_degC < bn(i+1));
    Ts_ch4_3a{i} = gs.ch4_eq1(gs.Soil_10_degC >= bn(i) & gs.Soil_10_degC < bn(i+1));
    Ts_ch4_3b{i} = gs.ch4_eq2(gs.Soil_10_degC >= bn(i) & gs.Soil_10_degC < bn(i+1));
    Ts_ghg_avg3a(i,1) = nanmean(Ts_ghg3a{i});
    Ts_ghg_avg3b(i,1) = nanmean(Ts_ghg3b{i});
    Ts_co2_avg3(i,1) = nanmean(Ts_co2_3{i});
    Ts_ch4_avg3a(i,1) = nanmean(Ts_ch4_3a{i});
    Ts_ch4_avg3b(i,1) = nanmean(Ts_ch4_3b{i});
    Ts_ghg_std3a(i,1) = nanstd(Ts_ghg3a{i});
    Ts_ghg_std3b(i,1) = nanstd(Ts_ghg3b{i});
    Ts_co2_std3(i,1) = nanstd(Ts_co2_3{i});
    Ts_ch4_std3a(i,1) = nanstd(Ts_ch4_3a{i});
    Ts_ch4_std3b(i,1) = nanstd(Ts_ch4_3b{i});
    Ts_ghg_se3a(i,1) = Ts_ghg_std3a(i,1)/sqrt(length(Ts_ghg3a{i}));
    Ts_ghg_se3b(i,1) = Ts_ghg_std3b(i,1)/sqrt(length(Ts_ghg3b{i}));
    Ts_co2_se3(i,1) = Ts_co2_std3(i,1)/sqrt(length(Ts_co2_3{i}));
    Ts_ch4_se3a(i,1) = Ts_ch4_std3a(i,1)/sqrt(length(Ts_ch4_3a{i}));
    Ts_ch4_se3b(i,1) = Ts_ch4_std3b(i,1)/sqrt(length(Ts_ch4_3b{i}));
end

bn = linspace(-0.5, 12, 16);
for i = 1:15
    bin_Ts4{i} = sf.Soil_10_degC(sf.Soil_10_degC >= bn(i) & sf.Soil_10_degC < bn(i+1));
    bin_Ts_avg4(i,1) = nanmean(bin_Ts4{i});
    Ts_ghg4a{i} = sf.ghg1(sf.Soil_10_degC >= bn(i) & sf.Soil_10_degC < bn(i+1));
    Ts_ghg4b{i} = sf.ghg2(sf.Soil_10_degC >= bn(i) & sf.Soil_10_degC < bn(i+1));
    Ts_co2_4{i} = sf.co2(sf.Soil_10_degC >= bn(i) & sf.Soil_10_degC < bn(i+1));
    Ts_ch4_4a{i} = sf.ch4_eq1(sf.Soil_10_degC >= bn(i) & sf.Soil_10_degC < bn(i+1));
    Ts_ch4_4b{i} = sf.ch4_eq2(sf.Soil_10_degC >= bn(i) & sf.Soil_10_degC < bn(i+1));
    Ts_ghg_avg4a(i,1) = nanmean(Ts_ghg4a{i});
    Ts_ghg_avg4b(i,1) = nanmean(Ts_ghg4b{i});
    Ts_co2_avg4(i,1) = nanmean(Ts_co2_4{i});
    Ts_ch4_avg4a(i,1) = nanmean(Ts_ch4_4a{i});
    Ts_ch4_avg4b(i,1) = nanmean(Ts_ch4_4b{i});
    Ts_ghg_std4a(i,1) = nanstd(Ts_ghg4a{i});
    Ts_ghg_std4b(i,1) = nanstd(Ts_ghg4b{i});
    Ts_co2_std4(i,1) = nanstd(Ts_co2_4{i});
    Ts_ch4_std4a(i,1) = nanstd(Ts_ch4_4a{i});
    Ts_ch4_std4b(i,1) = nanstd(Ts_ch4_4b{i});
    Ts_ghg_se4a(i,1) = Ts_ghg_std4a(i,1)/sqrt(length(Ts_ghg4a{i}));
    Ts_ghg_se4b(i,1) = Ts_ghg_std4b(i,1)/sqrt(length(Ts_ghg4b{i}));
    Ts_co2_se4(i,1) = Ts_co2_std4(i,1)/sqrt(length(Ts_co2_4{i}));
    Ts_ch4_se4a(i,1) = Ts_ch4_std4a(i,1)/sqrt(length(Ts_ch4_4a{i}));
    Ts_ch4_se4b(i,1) = Ts_ch4_std4b(i,1)/sqrt(length(Ts_ch4_4b{i}));
end

%%%%%%%%%%%%%%%%Tair responses%%%%%%%%%%%%%%%%%%
bn = linspace(-30, 15, 16);
for i = 1:15
    bin_Ta1{i} = w.Ta_degC(w.Ta_degC >= bn(i) & w.Ta_degC < bn(i+1));
    bin_Ta_avg1(i,1) = nanmean(bin_Ta1{i});
    Ta_ghg1a{i} = w.ghg1(w.Ta_degC >= bn(i) & w.Ta_degC < bn(i+1));
    Ta_ghg1b{i} = w.ghg2(w.Ta_degC >= bn(i) & w.Ta_degC < bn(i+1));
    Ta_co2_1{i} = w.co2(w.Ta_degC >= bn(i) & w.Ta_degC < bn(i+1));
    Ta_ch4_1a{i} = w.ch4_eq1(w.Ta_degC >= bn(i) & w.Ta_degC < bn(i+1));
    Ta_ch4_1b{i} = w.ch4_eq2(w.Ta_degC >= bn(i) & w.Ta_degC < bn(i+1));
    Ta_ghg_avg1a(i,1) = nanmean(Ta_ghg1a{i});
    Ta_ghg_avg1b(i,1) = nanmean(Ta_ghg1b{i});
    Ta_co2_avg1(i,1) = nanmean(Ta_co2_1{i}); 
    Ta_ch4_avg1a(i,1) = nanmean(Ta_ch4_1a{i});
    Ta_ch4_avg1b(i,1) = nanmean(Ta_ch4_1b{i});
    Ta_ghg_std1a(i,1) = nanstd(Ta_ghg1a{i});
    Ta_ghg_std1b(i,1) = nanstd(Ta_ghg1b{i});
    Ta_co2_std1(i,1) = nanstd(Ta_co2_1{i});
    Ta_ch4_std1a(i,1) = nanstd(Ta_ch4_1a{i});
    Ta_ch4_std1b(i,1) = nanstd(Ta_ch4_1b{i});
    Ta_ghg_se1a(i,1) = Ta_ghg_std1a(i,1)/sqrt(length(Ta_ghg1a{i}));
    Ta_ghg_se1b(i,1) = Ta_ghg_std1b(i,1)/sqrt(length(Ta_ghg1b{i}));
    Ta_co2_se1(i,1) = Ta_co2_std1(i,1)/sqrt(length(Ta_co2_1{i}));
    Ta_ch4_se1a(i,1) = Ta_ch4_std1a(i,1)/sqrt(length(Ta_ch4_1a{i}));
    Ta_ch4_se1b(i,1) = Ta_ch4_std1b(i,1)/sqrt(length(Ta_ch4_1b{i}));
end

bn = linspace(-15, 22, 16);
for i = 1:15
    bin_Ta2{i} = st.Ta_degC(st.Ta_degC >= bn(i) & st.Ta_degC < bn(i+1));
    bin_Ta_avg2(i,1) = nanmean(bin_Ta2{i});
    Ta_ghg2a{i} = st.ghg1(st.Ta_degC >= bn(i) & st.Ta_degC < bn(i+1));
    Ta_ghg2b{i} = st.ghg2(st.Ta_degC >= bn(i) & st.Ta_degC < bn(i+1));
    Ta_co2_2{i} = st.co2(st.Ta_degC >= bn(i) & st.Ta_degC < bn(i+1));
    Ta_ch4_2a{i} = st.ch4_eq1(st.Ta_degC >= bn(i) & st.Ta_degC < bn(i+1));
    Ta_ch4_2b{i} = st.ch4_eq2(st.Ta_degC >= bn(i) & st.Ta_degC < bn(i+1));
    Ta_ghg_avg2a(i,1) = nanmean(Ta_ghg2a{i});
    Ta_ghg_avg2b(i,1) = nanmean(Ta_ghg2b{i});
    Ta_co2_avg2(i,1) = nanmean(Ta_co2_2{i});
    Ta_ch4_avg2a(i,1) = nanmean(Ta_ch4_2a{i});
    Ta_ch4_avg2b(i,1) = nanmean(Ta_ch4_2b{i});
    Ta_ghg_std2a(i,1) = nanstd(Ta_ghg2a{i});
    Ta_ghg_std2b(i,1) = nanstd(Ta_ghg2b{i});
    Ta_co2_std2(i,1) = nanstd(Ta_co2_2{i});
    Ta_ch4_std2a(i,1) = nanstd(Ta_ch4_2a{i});
    Ta_ch4_std2b(i,1) = nanstd(Ta_ch4_2b{i});
    Ta_ghg_se2a(i,1) = Ta_ghg_std2a(i,1)/sqrt(length(Ta_ghg2a{i}));
    Ta_ghg_se2b(i,1) = Ta_ghg_std2b(i,1)/sqrt(length(Ta_ghg2b{i}));
    Ta_co2_se2(i,1) = Ta_co2_std2(i,1)/sqrt(length(Ta_co2_2{i}));
    Ta_ch4_se2a(i,1) = Ta_ch4_std2a(i,1)/sqrt(length(Ta_ch4_2a{i}));
    Ta_ch4_se2b(i,1) = Ta_ch4_std2b(i,1)/sqrt(length(Ta_ch4_2b{i}));
end

bn = linspace(-8, 25, 16);
for i = 1:15
    bin_Ta3{i} = gs.Ta_degC(gs.Ta_degC >= bn(i) & gs.Ta_degC < bn(i+1));
    bin_Ta_avg3(i,1) = nanmean(bin_Ta3{i});
    Ta_ghg3a{i} = gs.ghg1(gs.Ta_degC >= bn(i) & gs.Ta_degC < bn(i+1));
    Ta_ghg3b{i} = gs.ghg2(gs.Ta_degC >= bn(i) & gs.Ta_degC < bn(i+1));
    Ta_co2_3{i} = gs.co2(gs.Ta_degC >= bn(i) & gs.Ta_degC < bn(i+1));
    Ta_ch4_3a{i} = gs.ch4_eq1(gs.Ta_degC >= bn(i) & gs.Ta_degC < bn(i+1));
    Ta_ch4_3b{i} = gs.ch4_eq2(gs.Ta_degC >= bn(i) & gs.Ta_degC < bn(i+1));
    Ta_ghg_avg3a(i,1) = nanmean(Ta_ghg3a{i});
    Ta_ghg_avg3b(i,1) = nanmean(Ta_ghg3b{i});
    Ta_co2_avg3(i,1) = nanmean(Ta_co2_3{i});
    Ta_ch4_avg3a(i,1) = nanmean(Ta_ch4_3a{i});
    Ta_ch4_avg3b(i,1) = nanmean(Ta_ch4_3b{i});
    Ta_ghg_std3a(i,1) = nanstd(Ta_ghg3a{i});
    Ta_ghg_std3b(i,1) = nanstd(Ta_ghg3b{i});
    Ta_co2_std3(i,1) = nanstd(Ta_co2_3{i});
    Ta_ch4_std3a(i,1) = nanstd(Ta_ch4_3a{i});
    Ta_ch4_std3b(i,1) = nanstd(Ta_ch4_3b{i});
    Ta_ghg_se3a(i,1) = Ta_ghg_std3a(i,1)/sqrt(length(Ta_ghg3a{i}));
    Ta_ghg_se3b(i,1) = Ta_ghg_std3b(i,1)/sqrt(length(Ta_ghg3b{i}));
    Ta_co2_se3(i,1) = Ta_co2_std3(i,1)/sqrt(length(Ta_co2_3{i}));
    Ta_ch4_se3a(i,1) = Ta_ch4_std3a(i,1)/sqrt(length(Ta_ch4_3a{i}));
    Ta_ch4_se3b(i,1) = Ta_ch4_std3b(i,1)/sqrt(length(Ta_ch4_3b{i}));
end

bn = linspace(-20, 20, 16);
for i = 1:15
    bin_Ta4{i} = sf.Ta_degC(sf.Ta_degC >= bn(i) & sf.Ta_degC < bn(i+1));
    bin_Ta_avg4(i,1) = nanmean(bin_Ta4{i});
    Ta_ghg4a{i} = sf.ghg1(sf.Ta_degC >= bn(i) & sf.Ta_degC < bn(i+1));
    Ta_ghg4b{i} = sf.ghg2(sf.Ta_degC >= bn(i) & sf.Ta_degC < bn(i+1));
    Ta_co2_4{i} = sf.co2(sf.Ta_degC >= bn(i) & sf.Ta_degC < bn(i+1));
    Ta_ch4_4a{i} = sf.ch4_eq1(sf.Ta_degC >= bn(i) & sf.Ta_degC < bn(i+1));
    Ta_ch4_4b{i} = sf.ch4_eq2(sf.Ta_degC >= bn(i) & sf.Ta_degC < bn(i+1));
    Ta_ghg_avg4a(i,1) = nanmean(Ta_ghg4a{i});
    Ta_ghg_avg4b(i,1) = nanmean(Ta_ghg4b{i});
    Ta_co2_avg4(i,1) = nanmean(Ta_co2_4{i});
    Ta_ch4_avg4a(i,1) = nanmean(Ta_ch4_4a{i});
    Ta_ch4_avg4b(i,1) = nanmean(Ta_ch4_4b{i});
    Ta_ghg_std4a(i,1) = nanstd(Ta_ghg4a{i});
    Ta_ghg_std4b(i,1) = nanstd(Ta_ghg4b{i});
    Ta_co2_std4(i,1) = nanstd(Ta_co2_4{i});
    Ta_ch4_std4a(i,1) = nanstd(Ta_ch4_4a{i});
    Ta_ch4_std4b(i,1) = nanstd(Ta_ch4_4b{i});
    Ta_ghg_se4a(i,1) = Ta_ghg_std4a(i,1)/sqrt(length(Ta_ghg4a{i}));
    Ta_ghg_se4b(i,1) = Ta_ghg_std4b(i,1)/sqrt(length(Ta_ghg4b{i}));
    Ta_co2_se4(i,1) = Ta_co2_std4(i,1)/sqrt(length(Ta_co2_4{i}));
    Ta_ch4_se4a(i,1) = Ta_ch4_std4a(i,1)/sqrt(length(Ta_ch4_4a{i}));
    Ta_ch4_se4b(i,1) = Ta_ch4_std4b(i,1)/sqrt(length(Ta_ch4_4b{i}));
end

%%%%%%%%%%%%%%%VPD responses%%%%%%%%%%%%%%%%%%%%%
bn = linspace(0, 1.5, 16);
for i = 1:15
    bin_v1{i} = w.VPD_kPa(w.VPD_kPa >= bn(i) & w.VPD_kPa < bn(i+1));
    bin_v_avg1(i,1) = nanmean(bin_v1{i});
    v_ghg1a{i} = w.ghg1(w.VPD_kPa >= bn(i) & w.VPD_kPa < bn(i+1));
    v_ghg1b{i} = w.ghg2(w.VPD_kPa >= bn(i) & w.VPD_kPa < bn(i+1));
    v_co2_1{i} = w.co2(w.VPD_kPa >= bn(i) & w.VPD_kPa < bn(i+1));
    v_ch4_1a{i} = w.ch4_eq1(w.VPD_kPa >= bn(i) & w.VPD_kPa < bn(i+1));
    v_ch4_1b{i} = w.ch4_eq2(w.VPD_kPa >= bn(i) & w.VPD_kPa < bn(i+1));
    v_ghg_avg1a(i,1) = nanmean(v_ghg1a{i});
    v_ghg_avg1b(i,1) = nanmean(v_ghg1b{i});
    v_co2_avg1(i,1) = nanmean(v_co2_1{i}); 
    v_ch4_avg1a(i,1) = nanmean(v_ch4_1a{i});
    v_ch4_avg1b(i,1) = nanmean(v_ch4_1b{i});
    v_ghg_std1a(i,1) = nanstd(v_ghg1a{i});
    v_ghg_std1b(i,1) = nanstd(v_ghg1b{i});
    v_co2_std1(i,1) = nanstd(v_co2_1{i});
    v_ch4_std1a(i,1) = nanstd(v_ch4_1a{i});
    v_ch4_std1b(i,1) = nanstd(v_ch4_1b{i});
    v_ghg_se1a(i,1) = v_ghg_std1a(i,1)/sqrt(length(v_ghg1a{i}));
    v_ghg_se1b(i,1) = v_ghg_std1b(i,1)/sqrt(length(v_ghg1b{i}));
    v_co2_se1(i,1) = v_co2_std1(i,1)/sqrt(length(v_co2_1{i}));
    v_ch4_se1a(i,1) = v_ch4_std1a(i,1)/sqrt(length(v_ch4_1a{i}));
    v_ch4_se1b(i,1) = v_ch4_std1b(i,1)/sqrt(length(v_ch4_1b{i}));
end

bn = linspace(0, 2.3, 16);
for i = 1:15
    bin_v2{i} = st.VPD_kPa(st.VPD_kPa >= bn(i) & st.VPD_kPa < bn(i+1));
    bin_v_avg2(i,1) = nanmean(bin_v2{i});
    v_ghg2a{i} = st.ghg1(st.VPD_kPa >= bn(i) & st.VPD_kPa < bn(i+1));
    v_ghg2b{i} = st.ghg2(st.VPD_kPa >= bn(i) & st.VPD_kPa < bn(i+1));
    v_co2_2{i} = st.co2(st.VPD_kPa >= bn(i) & st.VPD_kPa < bn(i+1));
    v_ch4_2a{i} = st.ch4_eq1(st.VPD_kPa >= bn(i) & st.VPD_kPa < bn(i+1));
    v_ch4_2b{i} = st.ch4_eq2(st.VPD_kPa >= bn(i) & st.VPD_kPa < bn(i+1));
    v_ghg_avg2a(i,1) = nanmean(v_ghg2a{i});
    v_ghg_avg2b(i,1) = nanmean(v_ghg2b{i});
    v_co2_avg2(i,1) = nanmean(v_co2_2{i});
    v_ch4_avg2a(i,1) = nanmean(v_ch4_2a{i});
    v_ch4_avg2b(i,1) = nanmean(v_ch4_2b{i});
    v_ghg_std2a(i,1) = nanstd(v_ghg2a{i});
    v_ghg_std2b(i,1) = nanstd(v_ghg2b{i});
    v_co2_std2(i,1) = nanstd(v_co2_2{i});
    v_ch4_std2a(i,1) = nanstd(v_ch4_2a{i});
    v_ch4_std2b(i,1) = nanstd(v_ch4_2b{i});
    v_ghg_se2a(i,1) = v_ghg_std2a(i,1)/sqrt(length(v_ghg2a{i}));
    v_ghg_se2b(i,1) = v_ghg_std2b(i,1)/sqrt(length(v_ghg2b{i}));
    v_co2_se2(i,1) = v_co2_std2(i,1)/sqrt(length(v_co2_2{i}));
    v_ch4_se2a(i,1) = v_ch4_std2a(i,1)/sqrt(length(v_ch4_2a{i}));
    v_ch4_se2b(i,1) = v_ch4_std2b(i,1)/sqrt(length(v_ch4_2b{i}));
end

bn = linspace(0, 2.2, 16);
for i = 1:15
    bin_v3{i} = gs.VPD_kPa(gs.VPD_kPa >= bn(i) & gs.VPD_kPa < bn(i+1));
    bin_v_avg3(i,1) = nanmean(bin_v3{i});
    v_ghg3a{i} = gs.ghg1(gs.VPD_kPa >= bn(i) & gs.VPD_kPa < bn(i+1));
    v_ghg3b{i} = gs.ghg2(gs.VPD_kPa >= bn(i) & gs.VPD_kPa < bn(i+1));
    v_co2_3{i} = gs.co2(gs.VPD_kPa >= bn(i) & gs.VPD_kPa < bn(i+1));
    v_ch4_3a{i} = gs.ch4_eq1(gs.VPD_kPa >= bn(i) & gs.VPD_kPa < bn(i+1));
    v_ch4_3b{i} = gs.ch4_eq2(gs.VPD_kPa >= bn(i) & gs.VPD_kPa < bn(i+1));
    v_ghg_avg3a(i,1) = nanmean(v_ghg3a{i});
    v_ghg_avg3b(i,1) = nanmean(v_ghg3b{i});
    v_co2_avg3(i,1) = nanmean(v_co2_3{i});
    v_ch4_avg3a(i,1) = nanmean(v_ch4_3a{i});
    v_ch4_avg3b(i,1) = nanmean(v_ch4_3b{i});
    v_ghg_std3a(i,1) = nanstd(v_ghg3a{i});
    v_ghg_std3b(i,1) = nanstd(v_ghg3b{i});
    v_co2_std3(i,1) = nanstd(v_co2_3{i});
    v_ch4_std3a(i,1) = nanstd(v_ch4_3a{i});
    v_ch4_std3b(i,1) = nanstd(v_ch4_3b{i});
    v_ghg_se3a(i,1) = v_ghg_std3a(i,1)/sqrt(length(v_ghg3a{i}));
    v_ghg_se3b(i,1) = v_ghg_std3b(i,1)/sqrt(length(v_ghg3b{i}));
    v_co2_se3(i,1) = v_co2_std3(i,1)/sqrt(length(v_co2_3{i}));
    v_ch4_se3a(i,1) = v_ch4_std3a(i,1)/sqrt(length(v_ch4_3a{i}));
    v_ch4_se3b(i,1) = v_ch4_std3b(i,1)/sqrt(length(v_ch4_3b{i}));
end

bn = linspace(0, 1.8, 16);
for i = 1:15
    bin_v4{i} = sf.VPD_kPa(sf.VPD_kPa >= bn(i) & sf.VPD_kPa < bn(i+1));
    bin_v_avg4(i,1) = nanmean(bin_v4{i});
    v_ghg4a{i} = sf.ghg1(sf.VPD_kPa >= bn(i) & sf.VPD_kPa < bn(i+1));
    v_ghg4b{i} = sf.ghg2(sf.VPD_kPa >= bn(i) & sf.VPD_kPa < bn(i+1));
    v_co2_4{i} = sf.co2(sf.VPD_kPa >= bn(i) & sf.VPD_kPa < bn(i+1));
    v_ch4_4a{i} = sf.ch4_eq1(sf.VPD_kPa >= bn(i) & sf.VPD_kPa < bn(i+1));
    v_ch4_4b{i} = sf.ch4_eq2(sf.VPD_kPa >= bn(i) & sf.VPD_kPa < bn(i+1));
    v_ghg_avg4a(i,1) = nanmean(v_ghg4a{i});
    v_ghg_avg4b(i,1) = nanmean(v_ghg4b{i});
    v_co2_avg4(i,1) = nanmean(v_co2_4{i});
    v_ch4_avg4a(i,1) = nanmean(v_ch4_4a{i});
    v_ch4_avg4b(i,1) = nanmean(v_ch4_4b{i});
    v_ghg_std4a(i,1) = nanstd(v_ghg4a{i});
    v_ghg_std4b(i,1) = nanstd(v_ghg4b{i});
    v_co2_std4(i,1) = nanstd(v_co2_4{i});
    v_ch4_std4a(i,1) = nanstd(v_ch4_4a{i});
    v_ch4_std4b(i,1) = nanstd(v_ch4_4b{i});
    v_ghg_se4a(i,1) = v_ghg_std4a(i,1)/sqrt(length(v_ghg4a{i}));
    v_ghg_se4b(i,1) = v_ghg_std4b(i,1)/sqrt(length(v_ghg4b{i}));
    v_co2_se4(i,1) = v_co2_std4(i,1)/sqrt(length(v_co2_4{i}));
    v_ch4_se4a(i,1) = v_ch4_std4a(i,1)/sqrt(length(v_ch4_4a{i}));
    v_ch4_se4b(i,1) = v_ch4_std4b(i,1)/sqrt(length(v_ch4_4b{i}));
end

%%%%%%%%%%%%%%%%%%SWC responses%%%%%%%%%%%%%%%%%%%
bn = linspace(0.15, 0.4, 16);
for i = 1:15
    bin_sm1{i} = w.SWC_10_m3m3(w.SWC_10_m3m3 >= bn(i) & w.SWC_10_m3m3 < bn(i+1));
    bin_sm_avg1(i,1) = nanmean(bin_sm1{i});
    sm_ghg1a{i} = w.ghg1(w.SWC_10_m3m3 >= bn(i) & w.SWC_10_m3m3 < bn(i+1));
    sm_ghg1b{i} = w.ghg2(w.SWC_10_m3m3 >= bn(i) & w.SWC_10_m3m3 < bn(i+1));
    sm_co2_1{i} = w.co2(w.SWC_10_m3m3 >= bn(i) & w.SWC_10_m3m3 < bn(i+1));
    sm_ch4_1a{i} = w.ch4_eq1(w.SWC_10_m3m3 >= bn(i) & w.SWC_10_m3m3 < bn(i+1));
    sm_ch4_1b{i} = w.ch4_eq2(w.SWC_10_m3m3 >= bn(i) & w.SWC_10_m3m3 < bn(i+1));
    sm_ghg_avg1a(i,1) = nanmean(sm_ghg1a{i});
    sm_ghg_avg1b(i,1) = nanmean(sm_ghg1b{i});
    sm_co2_avg1(i,1) = nanmean(sm_co2_1{i}); 
    sm_ch4_avg1a(i,1) = nanmean(sm_ch4_1a{i});
    sm_ch4_avg1b(i,1) = nanmean(sm_ch4_1b{i});
    sm_ghg_std1a(i,1) = nanstd(sm_ghg1a{i});
    sm_ghg_std1b(i,1) = nanstd(sm_ghg1b{i});
    sm_co2_std1(i,1) = nanstd(sm_co2_1{i});
    sm_ch4_std1a(i,1) = nanstd(sm_ch4_1a{i});
    sm_ch4_std1b(i,1) = nanstd(sm_ch4_1b{i});
    sm_ghg_se1a(i,1) = sm_ghg_std1a(i,1)/sqrt(length(sm_ghg1a{i}));
    sm_ghg_se1b(i,1) = sm_ghg_std1b(i,1)/sqrt(length(sm_ghg1b{i}));
    sm_co2_se1(i,1) = sm_co2_std1(i,1)/sqrt(length(sm_co2_1{i}));
    sm_ch4_se1a(i,1) = sm_ch4_std1a(i,1)/sqrt(length(sm_ch4_1a{i}));
    sm_ch4_se1b(i,1) = sm_ch4_std1b(i,1)/sqrt(length(sm_ch4_1b{i}));
end

bn = linspace(0.2, 0.52, 16);
for i = 1:15
    bin_sm2{i} = st.SWC_10_m3m3(st.SWC_10_m3m3 >= bn(i) & st.SWC_10_m3m3 < bn(i+1));
    bin_sm_avg2(i,1) = nanmean(bin_sm2{i});
    sm_ghg2a{i} = st.ghg1(st.SWC_10_m3m3 >= bn(i) & st.SWC_10_m3m3 < bn(i+1));
    sm_ghg2b{i} = st.ghg2(st.SWC_10_m3m3 >= bn(i) & st.SWC_10_m3m3 < bn(i+1));
    sm_co2_2{i} = st.co2(st.SWC_10_m3m3 >= bn(i) & st.SWC_10_m3m3 < bn(i+1));
    sm_ch4_2a{i} = st.ch4_eq1(st.SWC_10_m3m3 >= bn(i) & st.SWC_10_m3m3 < bn(i+1));
    sm_ch4_2b{i} = st.ch4_eq2(st.SWC_10_m3m3 >= bn(i) & st.SWC_10_m3m3 < bn(i+1));
    sm_ghg_avg2a(i,1) = nanmean(sm_ghg2a{i});
    sm_ghg_avg2b(i,1) = nanmean(sm_ghg2b{i});
    sm_co2_avg2(i,1) = nanmean(sm_co2_2{i});
    sm_ch4_avg2a(i,1) = nanmean(sm_ch4_2a{i});
    sm_ch4_avg2b(i,1) = nanmean(sm_ch4_2b{i});
    sm_ghg_std2a(i,1) = nanstd(sm_ghg2a{i});
    sm_ghg_std2b(i,1) = nanstd(sm_ghg2b{i});
    sm_co2_std2(i,1) = nanstd(sm_co2_2{i});
    sm_ch4_std2a(i,1) = nanstd(sm_ch4_2a{i});
    sm_ch4_std2b(i,1) = nanstd(sm_ch4_2b{i});
    sm_ghg_se2a(i,1) = sm_ghg_std2a(i,1)/sqrt(length(sm_ghg2a{i}));
    sm_ghg_se2b(i,1) = sm_ghg_std2b(i,1)/sqrt(length(sm_ghg2b{i}));
    sm_co2_se2(i,1) = sm_co2_std2(i,1)/sqrt(length(sm_co2_2{i}));
    sm_ch4_se2a(i,1) = sm_ch4_std2a(i,1)/sqrt(length(sm_ch4_2a{i}));
    sm_ch4_se2b(i,1) = sm_ch4_std2b(i,1)/sqrt(length(sm_ch4_2b{i}));
end

bn = linspace(0.36, 0.54, 16);
for i = 1:15
    bin_sm3{i} = gs.SWC_10_m3m3(gs.SWC_10_m3m3 >= bn(i) & gs.SWC_10_m3m3 < bn(i+1));
    bin_sm_avg3(i,1) = nanmean(bin_sm3{i});
    sm_ghg3a{i} = gs.ghg1(gs.SWC_10_m3m3 >= bn(i) & gs.SWC_10_m3m3 < bn(i+1));
    sm_ghg3b{i} = gs.ghg2(gs.SWC_10_m3m3 >= bn(i) & gs.SWC_10_m3m3 < bn(i+1));
    sm_co2_3{i} = gs.co2(gs.SWC_10_m3m3 >= bn(i) & gs.SWC_10_m3m3 < bn(i+1));
    sm_ch4_3a{i} = gs.ch4_eq1(gs.SWC_10_m3m3 >= bn(i) & gs.SWC_10_m3m3 < bn(i+1));
    sm_ch4_3b{i} = gs.ch4_eq2(gs.SWC_10_m3m3 >= bn(i) & gs.SWC_10_m3m3 < bn(i+1));
    sm_ghg_avg3a(i,1) = nanmean(sm_ghg3a{i});
    sm_ghg_avg3b(i,1) = nanmean(sm_ghg3b{i});
    sm_co2_avg3(i,1) = nanmean(sm_co2_3{i});
    sm_ch4_avg3a(i,1) = nanmean(sm_ch4_3a{i});
    sm_ch4_avg3b(i,1) = nanmean(sm_ch4_3b{i});
    sm_ghg_std3a(i,1) = nanstd(sm_ghg3a{i});
    sm_ghg_std3b(i,1) = nanstd(sm_ghg3b{i});
    sm_co2_std3(i,1) = nanstd(sm_co2_3{i});
    sm_ch4_std3a(i,1) = nanstd(sm_ch4_3a{i});
    sm_ch4_std3b(i,1) = nanstd(sm_ch4_3b{i});
    sm_ghg_se3a(i,1) = sm_ghg_std3a(i,1)/sqrt(length(sm_ghg3a{i}));
    sm_ghg_se3b(i,1) = sm_ghg_std3b(i,1)/sqrt(length(sm_ghg3b{i}));
    sm_co2_se3(i,1) = sm_co2_std3(i,1)/sqrt(length(sm_co2_3{i}));
    sm_ch4_se3a(i,1) = sm_ch4_std3a(i,1)/sqrt(length(sm_ch4_3a{i}));
    sm_ch4_se3b(i,1) = sm_ch4_std3b(i,1)/sqrt(length(sm_ch4_3b{i}));
end

bn = linspace(0.29, 0.54, 16);
for i = 1:15
    bin_sm4{i} = sf.SWC_10_m3m3(sf.SWC_10_m3m3 >= bn(i) & sf.SWC_10_m3m3 < bn(i+1));
    bin_sm_avg4(i,1) = nanmean(bin_sm4{i});
    sm_ghg4a{i} = sf.ghg1(sf.SWC_10_m3m3 >= bn(i) & sf.SWC_10_m3m3 < bn(i+1));
    sm_ghg4b{i} = sf.ghg2(sf.SWC_10_m3m3 >= bn(i) & sf.SWC_10_m3m3 < bn(i+1));
    sm_co2_4{i} = sf.co2(sf.SWC_10_m3m3 >= bn(i) & sf.SWC_10_m3m3 < bn(i+1));
    sm_ch4_4a{i} = sf.ch4_eq1(sf.SWC_10_m3m3 >= bn(i) & sf.SWC_10_m3m3 < bn(i+1));
    sm_ch4_4b{i} = sf.ch4_eq2(sf.SWC_10_m3m3 >= bn(i) & sf.SWC_10_m3m3 < bn(i+1));
    sm_ghg_avg4a(i,1) = nanmean(sm_ghg4a{i});
    sm_ghg_avg4b(i,1) = nanmean(sm_ghg4b{i});
    sm_co2_avg4(i,1) = nanmean(sm_co2_4{i});
    sm_ch4_avg4a(i,1) = nanmean(sm_ch4_4a{i});
    sm_ch4_avg4b(i,1) = nanmean(sm_ch4_4b{i});
    sm_ghg_std4a(i,1) = nanstd(sm_ghg4a{i});
    sm_ghg_std4b(i,1) = nanstd(sm_ghg4b{i});
    sm_co2_std4(i,1) = nanstd(sm_co2_4{i});
    sm_ch4_std4a(i,1) = nanstd(sm_ch4_4a{i});
    sm_ch4_std4b(i,1) = nanstd(sm_ch4_4b{i});
    sm_ghg_se4a(i,1) = sm_ghg_std4a(i,1)/sqrt(length(sm_ghg4a{i}));
    sm_ghg_se4b(i,1) = sm_ghg_std4b(i,1)/sqrt(length(sm_ghg4b{i}));
    sm_co2_se4(i,1) = sm_co2_std4(i,1)/sqrt(length(sm_co2_4{i}));
    sm_ch4_se4a(i,1) = sm_ch4_std4a(i,1)/sqrt(length(sm_ch4_4a{i}));
    sm_ch4_se4b(i,1) = sm_ch4_std4b(i,1)/sqrt(length(sm_ch4_4b{i}));
end

%%%%%%%%%%%%%%PPT response%%%%%%%%%%%%%%%
w_sum = retime(w, 'daily', 'sum');
w_avg = retime(w, 'daily', 'mean');

bn = linspace(0, 4.5, 11);
for i = 1:10
    bin_p1{i} = w_sum.PPT_mm(w_sum.PPT_mm >= bn(i) & w_sum.PPT_mm < bn(i+1));
    bin_p_avg1(i,1) = nanmean(bin_p1{i});
    p_ghg1a{i} = w_avg.ghg1(w_sum.PPT_mm >= bn(i) & w_sum.PPT_mm < bn(i+1));
    p_ghg1b{i} = w_avg.ghg2(w_sum.PPT_mm >= bn(i) & w_sum.PPT_mm < bn(i+1));
    p_co2_1{i} = w_avg.co2(w_sum.PPT_mm >= bn(i) & w_sum.PPT_mm < bn(i+1));
    p_ch4_1a{i} = w_avg.ch4_eq1(w_sum.PPT_mm >= bn(i) & w_sum.PPT_mm < bn(i+1));
    p_ch4_1b{i} = w_avg.ch4_eq2(w_sum.PPT_mm >= bn(i) & w_sum.PPT_mm < bn(i+1));
    p_ghg_avg1a(i,1) = nanmean(p_ghg1a{i});
    p_ghg_avg1b(i,1) = nanmean(p_ghg1b{i});
    p_co2_avg1(i,1) = nanmean(p_co2_1{i}); 
    p_ch4_avg1a(i,1) = nanmean(p_ch4_1a{i});
    p_ch4_avg1b(i,1) = nanmean(p_ch4_1b{i});
    p_ghg_std1a(i,1) = nanstd(p_ghg1a{i});
    p_ghg_std1b(i,1) = nanstd(p_ghg1b{i});
    p_co2_std1(i,1) = nanstd(p_co2_1{i});
    p_ch4_std1a(i,1) = nanstd(p_ch4_1a{i});
    p_ch4_std1b(i,1) = nanstd(p_ch4_1b{i});
    p_ghg_se1a(i,1) = p_ghg_std1a(i,1)/sqrt(length(p_ghg1a{i}));
    p_ghg_se1b(i,1) = p_ghg_std1b(i,1)/sqrt(length(p_ghg1b{i}));
    p_co2_se1(i,1) = p_co2_std1(i,1)/sqrt(length(p_co2_1{i}));
    p_ch4_se1a(i,1) = p_ch4_std1a(i,1)/sqrt(length(p_ch4_1a{i}));
    p_ch4_se1b(i,1) = p_ch4_std1b(i,1)/sqrt(length(p_ch4_1b{i}));
end


st_sum = retime(st, 'daily', 'sum');
st_avg = retime(st, 'daily', 'mean');

bn = linspace(0, 11, 11);
for i = 1:10
    bin_p2{i} = st_sum.PPT_mm(st_sum.PPT_mm >= bn(i) & st_sum.PPT_mm < bn(i+1));
    bin_p_avg2(i,1) = nanmean(bin_p2{i});
    p_ghg2a{i} = st_avg.ghg1(st_sum.PPT_mm >= bn(i) & st_sum.PPT_mm < bn(i+1));
    p_ghg2b{i} = st_avg.ghg2(st_sum.PPT_mm >= bn(i) & st_sum.PPT_mm < bn(i+1));
    p_co2_2{i} = st_avg.co2(st_sum.PPT_mm >= bn(i) & st_sum.PPT_mm < bn(i+1));
    p_ch4_2a{i} = st_avg.ch4_eq1(st_sum.PPT_mm >= bn(i) & st_sum.PPT_mm < bn(i+1));
    p_ch4_2b{i} = st_avg.ch4_eq2(st_sum.PPT_mm >= bn(i) & st_sum.PPT_mm < bn(i+1));
    p_ghg_avg2a(i,1) = nanmean(p_ghg2a{i});
    p_ghg_avg2b(i,1) = nanmean(p_ghg2b{i});
    p_co2_avg2(i,1) = nanmean(p_co2_2{i});
    p_ch4_avg2a(i,1) = nanmean(p_ch4_2a{i});
    p_ch4_avg2b(i,1) = nanmean(p_ch4_2b{i});
    p_ghg_std2a(i,1) = nanstd(p_ghg2a{i});
    p_ghg_std2b(i,1) = nanstd(p_ghg2b{i});
    p_co2_std2(i,1) = nanstd(p_co2_2{i});
    p_ch4_std2a(i,1) = nanstd(p_ch4_2a{i});
    p_ch4_std2b(i,1) = nanstd(p_ch4_2b{i});
    p_ghg_se2a(i,1) = p_ghg_std2a(i,1)/sqrt(length(p_ghg2a{i}));
    p_ghg_se2b(i,1) = p_ghg_std2b(i,1)/sqrt(length(p_ghg2b{i}));
    p_co2_se2(i,1) = p_co2_std2(i,1)/sqrt(length(p_co2_2{i}));
    p_ch4_se2a(i,1) = p_ch4_std2a(i,1)/sqrt(length(p_ch4_2a{i}));
    p_ch4_se2b(i,1) = p_ch4_std2b(i,1)/sqrt(length(p_ch4_2b{i}));
end


gs_sum = retime(gs, 'daily', 'sum');
gs_avg = retime(gs, 'daily', 'mean');

bn = linspace(0, 35, 11);
for i = 1:10
    bin_p3{i} = gs_sum.PPT_mm(gs_sum.PPT_mm >= bn(i) & gs_sum.PPT_mm < bn(i+1));
    bin_p_avg3(i,1) = nanmean(bin_p3{i});
    p_ghg3a{i} = gs_avg.ghg1(gs_sum.PPT_mm >= bn(i) & gs_sum.PPT_mm < bn(i+1));
    p_ghg3b{i} = gs_avg.ghg2(gs_sum.PPT_mm >= bn(i) & gs_sum.PPT_mm < bn(i+1));
    p_co2_3{i} = gs_avg.co2(gs_sum.PPT_mm >= bn(i) & gs_sum.PPT_mm < bn(i+1));
    p_ch4_3a{i} = gs_avg.ch4_eq1(gs_sum.PPT_mm >= bn(i) & gs_sum.PPT_mm < bn(i+1));
    p_ch4_3b{i} = gs_avg.ch4_eq2(gs_sum.PPT_mm >= bn(i) & gs_sum.PPT_mm < bn(i+1));
    p_ghg_avg3a(i,1) = nanmean(p_ghg3a{i});
    p_ghg_avg3b(i,1) = nanmean(p_ghg3b{i});
    p_co2_avg3(i,1) = nanmean(p_co2_3{i});
    p_ch4_avg3a(i,1) = nanmean(p_ch4_3a{i});
    p_ch4_avg3b(i,1) = nanmean(p_ch4_3b{i});
    p_ghg_std3a(i,1) = nanstd(p_ghg3a{i});
    p_ghg_std3b(i,1) = nanstd(p_ghg3b{i});
    p_co2_std3(i,1) = nanstd(p_co2_3{i});
    p_ch4_std3a(i,1) = nanstd(p_ch4_3a{i});
    p_ch4_std3b(i,1) = nanstd(p_ch4_3b{i});
    p_ghg_se3a(i,1) = p_ghg_std3a(i,1)/sqrt(length(p_ghg3a{i}));
    p_ghg_se3b(i,1) = p_ghg_std3b(i,1)/sqrt(length(p_ghg3b{i}));
    p_co2_se3(i,1) = p_co2_std3(i,1)/sqrt(length(p_co2_3{i}));
    p_ch4_se3a(i,1) = p_ch4_std3a(i,1)/sqrt(length(p_ch4_3a{i}));
    p_ch4_se3b(i,1) = p_ch4_std3b(i,1)/sqrt(length(p_ch4_3b{i}));
end


sf_sum = retime(sf, 'daily', 'sum');
sf_avg = retime(sf, 'daily', 'mean');

bn = linspace(0, 8, 11);
for i = 1:10
    bin_p4{i} = sf.PPT_mm(sf.PPT_mm >= bn(i) & sf.PPT_mm < bn(i+1));
    bin_p_avg4(i,1) = nanmean(bin_p4{i});
    p_ghg4a{i} = sf.ghg1(sf.PPT_mm >= bn(i) & sf.PPT_mm < bn(i+1));
    p_ghg4b{i} = sf.ghg2(sf.PPT_mm >= bn(i) & sf.PPT_mm < bn(i+1));
    p_co2_4{i} = sf.co2(sf.PPT_mm >= bn(i) & sf.PPT_mm < bn(i+1));
    p_ch4_4a{i} = sf.ch4_eq1(sf.PPT_mm >= bn(i) & sf.PPT_mm < bn(i+1));
    p_ch4_4b{i} = sf.ch4_eq2(sf.PPT_mm >= bn(i) & sf.PPT_mm < bn(i+1));
    p_ghg_avg4a(i,1) = nanmean(p_ghg4a{i});
    p_ghg_avg4b(i,1) = nanmean(p_ghg4b{i});
    p_co2_avg4(i,1) = nanmean(p_co2_4{i});
    p_ch4_avg4a(i,1) = nanmean(p_ch4_4a{i});
    p_ch4_avg4b(i,1) = nanmean(p_ch4_4b{i});
    p_ghg_std4a(i,1) = nanstd(p_ghg4a{i});
    p_ghg_std4b(i,1) = nanstd(p_ghg4b{i});
    p_co2_std4(i,1) = nanstd(p_co2_4{i});
    p_ch4_std4a(i,1) = nanstd(p_ch4_4a{i});
    p_ch4_std4b(i,1) = nanstd(p_ch4_4b{i});
    p_ghg_se4a(i,1) = p_ghg_std4a(i,1)/sqrt(length(p_ghg4a{i}));
    p_ghg_se4b(i,1) = p_ghg_std4b(i,1)/sqrt(length(p_ghg4b{i}));
    p_co2_se4(i,1) = p_co2_std4(i,1)/sqrt(length(p_co2_4{i}));
    p_ch4_se4a(i,1) = p_ch4_std4a(i,1)/sqrt(length(p_ch4_4a{i}));
    p_ch4_se4b(i,1) = p_ch4_std4b(i,1)/sqrt(length(p_ch4_4b{i}));
end

figure(8)
set(gcf, 'Position', [100 100 1200 650]);
subplot(2,3,1);
hold on; box on;
errorbar(bin_Rg_avg1, Rg_ghg_avg1a, Rg_ghg_se1a, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_Rg_avg2, Rg_ghg_avg2a, Rg_ghg_se2a, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_Rg_avg3, Rg_ghg_avg3a, Rg_ghg_se3a, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_Rg_avg4, Rg_ghg_avg4a, Rg_ghg_se4a, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
ylim([-0.4 0.4]);
xlim([0 1500]);
xlabel('R_g (W m^{-2})');
ylabel('net CO_2-eq flux (mg CO_2-eq m^{-2} s^{-1})');
legend('W', 'ST', 'GS', 'SF');
title('(a)');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);

subplot(2,3,2);
hold on; box on;
errorbar(bin_Ta_avg1, Ta_ghg_avg1a, Ta_ghg_se1a, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_Ta_avg2, Ta_ghg_avg2a, Ta_ghg_se2a, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_Ta_avg3, Ta_ghg_avg3a, Ta_ghg_se3a, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_Ta_avg4, Ta_ghg_avg4a, Ta_ghg_se4a, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
ylim([-0.3 0.3]);
xlim([-30 30]);
xlabel('T_{air} (\circC)');
ylabel('net CO_2-eq flux (mg CO_2-eq m^{-2} s^{-1})');
title('(b)');
set(gca, 'YTick', -0.3:0.1:0.3, 'XTick', -30:10:30, 'FontName', 'Times New Roman', 'FontSize', 12);

subplot(2,3,3);
hold on; box on;
errorbar(bin_p_avg1(~isnan(p_ghg_avg1a)), p_ghg_avg1a(~isnan(p_ghg_avg1a)), p_ghg_se1a(~isnan(p_ghg_avg1a)),...
    'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_p_avg2(~isnan(p_ghg_avg2a)), p_ghg_avg2a(~isnan(p_ghg_avg2a)), p_ghg_se2a(~isnan(p_ghg_avg2a)),...
    'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_p_avg3(~isnan(p_ghg_avg3a)), p_ghg_avg3a(~isnan(p_ghg_avg3a)), p_ghg_se3a(~isnan(p_ghg_avg3a)),...
    'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_p_avg4(~isnan(p_ghg_avg4a)), p_ghg_avg4a(~isnan(p_ghg_avg4a)), p_ghg_se4a(~isnan(p_ghg_avg4a)),...
    'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
ylim([-0.2 0.2]);
xlabel('PPT (mm)');
ylabel('net CO_2-eq flux (mg CO_2-eq m^{-2} s^{-1})');
title('(c)');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);

subplot(2,3,4);
hold on; box on;
errorbar(bin_v_avg1, v_ghg_avg1a, v_ghg_se1a, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_v_avg2, v_ghg_avg2a, v_ghg_se2a, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_v_avg3, v_ghg_avg3a, v_ghg_se3a, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_v_avg4, v_ghg_avg4a, v_ghg_se4a, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
ylim([-0.3 0.3]);
xlim([0 2.5]);
xlabel('VPD (kPa)');
ylabel('net CO_2-eq flux (mg CO_2-eq m^{-2} s^{-1})');
title('(d)');
set(gca, 'YTick', -0.3:0.1:0.3, 'XTick', 0:0.5:2.5, 'FontName', 'Times New Roman', 'FontSize', 12);

subplot(2,3,5);
hold on; box on;
errorbar(bin_Ts_avg1, Ts_ghg_avg1a, Ts_ghg_se1a, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_Ts_avg2, Ts_ghg_avg2a, Ts_ghg_se2a, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_Ts_avg3, Ts_ghg_avg3a, Ts_ghg_se3a, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_Ts_avg4, Ts_ghg_avg4a, Ts_ghg_se4a, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
ylim([-0.2 0.2]);
xlim([-10 30]);
xlabel('T_{soil} 10 cm (\circC)');
ylabel('net CO_2-eq flux (mg CO_2-eq m^{-2} s^{-1})');
title('(e)');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);

subplot(2,3,6);
hold on; box on;
errorbar(bin_sm_avg1, sm_ghg_avg1a, sm_ghg_se1a, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_sm_avg2, sm_ghg_avg2a, sm_ghg_se2a, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_sm_avg3, sm_ghg_avg3a, sm_ghg_se3a, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_sm_avg4, sm_ghg_avg4a, sm_ghg_se4a, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
ylim([-0.2 0.2]);
xlim([0.1 0.6]);
xlabel('SWC 10 cm (m^{3} m^{-3})');
ylabel('net CO_2-eq flux (mg CO_2-eq m^{-2} s^{-1})');
title('(f)');
set(gca, 'XTick', 0.1:0.1:0.6, 'FontName', 'Times New Roman', 'FontSize', 12);


figure(9)
set(gcf, 'Position', [100 100 1200 650]);
subplot(2,3,1);
hold on; box on;
errorbar(bin_Rg_avg1, Rg_ghg_avg1b, Rg_ghg_se1b, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_Rg_avg2, Rg_ghg_avg2b, Rg_ghg_se2b, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_Rg_avg3, Rg_ghg_avg3b, Rg_ghg_se3b, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_Rg_avg4, Rg_ghg_avg4b, Rg_ghg_se4b, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
ylim([-0.4 0.4]);
xlim([0 1500]);
xlabel('R_g (W m^{-2})');
ylabel('net CO_2-eq flux (mg CO_2-eq m^{-2} s^{-1})');
legend('W', 'ST', 'GS', 'SF');
title('(a)');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);

subplot(2,3,2);
hold on; box on;
errorbar(bin_Ta_avg1, Ta_ghg_avg1b, Ta_ghg_se1b, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_Ta_avg2, Ta_ghg_avg2b, Ta_ghg_se2b, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_Ta_avg3, Ta_ghg_avg3b, Ta_ghg_se3b, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_Ta_avg4, Ta_ghg_avg4b, Ta_ghg_se4b, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
ylim([-0.3 0.3]);
xlim([-30 30]);
xlabel('T_{air} (\circC)');
ylabel('net CO_2-eq flux (mg CO_2-eq m^{-2} s^{-1})');
title('(b)');
set(gca, 'YTick', -0.3:0.1:0.3, 'XTick', -30:10:30, 'FontName', 'Times New Roman', 'FontSize', 12);

subplot(2,3,3);
hold on; box on;
errorbar(bin_p_avg1(~isnan(p_ghg_avg1b)), p_ghg_avg1b(~isnan(p_ghg_avg1b)), p_ghg_se1b(~isnan(p_ghg_avg1b)),...
    'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_p_avg2(~isnan(p_ghg_avg2b)), p_ghg_avg2b(~isnan(p_ghg_avg2b)), p_ghg_se2b(~isnan(p_ghg_avg2b)),...
    'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_p_avg3(~isnan(p_ghg_avg3b)), p_ghg_avg3b(~isnan(p_ghg_avg3b)), p_ghg_se3b(~isnan(p_ghg_avg3b)),...
    'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_p_avg4(~isnan(p_ghg_avg4b)), p_ghg_avg4b(~isnan(p_ghg_avg4b)), p_ghg_se4b(~isnan(p_ghg_avg4b)),...
    'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
ylim([-0.2 0.2]);
xlabel('PPT (mm)');
ylabel('net CO_2-eq flux (mg CO_2-eq m^{-2} s^{-1})');
title('(c)');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);

subplot(2,3,4);
hold on; box on;
errorbar(bin_v_avg1, v_ghg_avg1b, v_ghg_se1b, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_v_avg2, v_ghg_avg2b, v_ghg_se2b, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_v_avg3, v_ghg_avg3b, v_ghg_se3b, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_v_avg4, v_ghg_avg4b, v_ghg_se4b, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
ylim([-0.3 0.3]);
xlim([0 2.5]);
xlabel('VPD (kPa)');
ylabel('net CO_2-eq flux (mg CO_2-eq m^{-2} s^{-1})');
title('(d)');
set(gca, 'YTick', -0.3:0.1:0.3, 'XTick', 0:0.5:2.5, 'FontName', 'Times New Roman', 'FontSize', 12);

subplot(2,3,5);
hold on; box on;
errorbar(bin_Ts_avg1, Ts_ghg_avg1b, Ts_ghg_se1b, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_Ts_avg2, Ts_ghg_avg2b, Ts_ghg_se2b, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_Ts_avg3, Ts_ghg_avg3b, Ts_ghg_se3b, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_Ts_avg4, Ts_ghg_avg4b, Ts_ghg_se4b, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
ylim([-0.1 0.2]);
xlim([-10 30]);
xlabel('T_{soil} 10 cm (\circC)');
ylabel('net CO_2-eq flux (mg CO_2-eq m^{-2} s^{-1})');
title('(e)');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);

subplot(2,3,6);
hold on; box on;
errorbar(bin_sm_avg1, sm_ghg_avg1b, sm_ghg_se1b, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_sm_avg2, sm_ghg_avg2b, sm_ghg_se2b, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_sm_avg3, sm_ghg_avg3b, sm_ghg_se3b, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_sm_avg4, sm_ghg_avg4b, sm_ghg_se4b, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
ylim([-0.2 0.2]);
xlim([0.1 0.6]);
xlabel('SWC 10 cm (m^{3} m^{-3})');
ylabel('net CO_2-eq flux (mg CO_2-eq m^{-2} s^{-1})');
title('(f)');
set(gca, 'XTick', 0.1:0.1:0.6, 'FontName', 'Times New Roman', 'FontSize', 12);


figure(10)
set(gcf, 'Position', [100 100 1200 650]);
subplot(2,3,1);
hold on; box on;
errorbar(bin_Rg_avg1, Rg_co2_avg1, Rg_co2_se1, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_Rg_avg2, Rg_co2_avg2, Rg_co2_se2, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_Rg_avg3, Rg_co2_avg3, Rg_co2_se3, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_Rg_avg4, Rg_co2_avg4, Rg_co2_se4, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
ylim([-0.4 0.4]);
xlim([0 1500]);
xlabel('R_g (W m^{-2})');
ylabel('CO_2 flux (mg CO_2 m^{-2} s^{-1})');
legend('W', 'ST', 'GS', 'SF');
title('(a)');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);

subplot(2,3,2);
hold on; box on;
errorbar(bin_Ta_avg1, Ta_co2_avg1, Ta_co2_se1, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_Ta_avg2, Ta_co2_avg2, Ta_co2_se2, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_Ta_avg3, Ta_co2_avg3, Ta_co2_se3, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_Ta_avg4, Ta_co2_avg4, Ta_co2_se4, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
ylim([-0.4 0.4]);
xlim([-30 30]);
xlabel('T_{air} (\circC)');
ylabel('CO_2 flux (mg CO_2 m^{-2} s^{-1})');
title('(b)');
set(gca, 'YTick', -0.4:0.2:0.4, 'XTick', -30:10:30, 'FontName', 'Times New Roman', 'FontSize', 12);

subplot(2,3,3);
hold on; box on;
errorbar(bin_p_avg1(~isnan(p_co2_avg1)), p_co2_avg1(~isnan(p_co2_avg1)), p_co2_se1(~isnan(p_co2_avg1)),...
    'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_p_avg2(~isnan(p_co2_avg2)), p_co2_avg2(~isnan(p_co2_avg2)), p_co2_se2(~isnan(p_co2_avg2)),...
    'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_p_avg3(~isnan(p_co2_avg3)), p_co2_avg3(~isnan(p_co2_avg3)), p_co2_se3(~isnan(p_co2_avg3)),...
    'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_p_avg4(~isnan(p_co2_avg4)), p_co2_avg4(~isnan(p_co2_avg4)), p_co2_se4(~isnan(p_co2_avg4)),...
    'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
ylim([-0.2 0.2]);
xlabel('PPT (mm)');
ylabel('CO_2 flux (mg CO_2 m^{-2} s^{-1})');
title('(c)');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);

subplot(2,3,4);
hold on; box on;
errorbar(bin_v_avg1, v_co2_avg1, v_co2_se1, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_v_avg2, v_co2_avg2, v_co2_se2, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_v_avg3, v_co2_avg3, v_co2_se3, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_v_avg4, v_co2_avg4, v_co2_se4, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
ylim([-0.4 0.2]);
xlim([0 2.5]);
xlabel('VPD (kPa)');
ylabel('CO_2 flux (mg CO_2 m^{-2} s^{-1})');
title('(d)');
set(gca, 'XTick', 0:0.5:2.5, 'FontName', 'Times New Roman', 'FontSize', 12);

subplot(2,3,5);
hold on; box on;
errorbar(bin_Ts_avg1, Ts_co2_avg1, Ts_co2_se1, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_Ts_avg2, Ts_co2_avg2, Ts_co2_se2, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_Ts_avg3, Ts_co2_avg3, Ts_co2_se3, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_Ts_avg4, Ts_co2_avg4, Ts_co2_se4, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
ylim([-0.3 0.1]);
xlim([-10 30]);
xlabel('T_{soil} 10 cm (\circC)');
ylabel('CO_2 flux (mg CO_2 m^{-2} s^{-1})');
title('(e)');
set(gca, 'YTick', -0.3:0.1:0.1, 'FontName', 'Times New Roman', 'FontSize', 12);

subplot(2,3,6);
hold on; box on;
errorbar(bin_sm_avg1, sm_co2_avg1, sm_co2_se1, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_sm_avg2, sm_co2_avg2, sm_co2_se2, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_sm_avg3, sm_co2_avg3, sm_co2_se3, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_sm_avg4, sm_co2_avg4, sm_co2_se4, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
ylim([-0.2 0.1]);
xlim([0.1 0.6]);
xlabel('SWC 10 cm (m^{3} m^{-3})');
ylabel('CO_2 flux (mg CO_2 m^{-2} s^{-1})');
title('(f)');
set(gca, 'XTick', 0.1:0.1:0.6, 'FontName', 'Times New Roman', 'FontSize', 12);


figure(11)
set(gcf, 'Position', [100 100 1200 650]);
subplot(2,3,1);
hold on; box on;
errorbar(bin_Rg_avg1, Rg_ch4_avg1a, Rg_ch4_se1a, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_Rg_avg2, Rg_ch4_avg2a, Rg_ch4_se2a, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_Rg_avg3, Rg_ch4_avg3a, Rg_ch4_se3a, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_Rg_avg4, Rg_ch4_avg4a, Rg_ch4_se4a, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
ylim([0 0.1]);
xlim([0 1500]);
xlabel('R_g (W m^{-2})');
ylabel('CH_4 flux (mg CO_2-eq m^{-2} s^{-1})');
legend('W', 'ST', 'GS', 'SF', 'Location', 'northeast', 'Orientation', 'horizontal', 'FontSize' ,10, 'NumColumns', 2);
title('(a)');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);

subplot(2,3,2);
hold on; box on;
errorbar(bin_Ta_avg1, Ta_ch4_avg1a, Ta_ch4_se1a, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_Ta_avg2, Ta_ch4_avg2a, Ta_ch4_se2a, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_Ta_avg3, Ta_ch4_avg3a, Ta_ch4_se3a, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_Ta_avg4, Ta_ch4_avg4a, Ta_ch4_se4a, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
ylim([-0.05 0.15]);
xlim([-30 30]);
xlabel('T_{air} (\circC)');
ylabel('CH_4 flux (mg CO_2-eq m^{-2} s^{-1})');
title('(b)');
set(gca, 'XTick', -30:10:30, 'FontName', 'Times New Roman', 'FontSize', 12);

subplot(2,3,3);
hold on; box on;
errorbar(bin_p_avg1(~isnan(p_ch4_avg1a)), p_ch4_avg1a(~isnan(p_ch4_avg1a)), p_ch4_se1a(~isnan(p_ch4_avg1a)),...
    'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_p_avg2(~isnan(p_ch4_avg2a)), p_ch4_avg2a(~isnan(p_ch4_avg2a)), p_ch4_se2a(~isnan(p_ch4_avg2a)),...
    'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_p_avg3(~isnan(p_ch4_avg3a)), p_ch4_avg3a(~isnan(p_ch4_avg3a)), p_ch4_se3a(~isnan(p_ch4_avg3a)),...
    'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_p_avg4(~isnan(p_ch4_avg4a)), p_ch4_avg4a(~isnan(p_ch4_avg4a)), p_ch4_se4a(~isnan(p_ch4_avg4a)),...
    'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
ylim([0 0.15]);
xlabel('PPT (mm)');
ylabel('CH_4 flux (mg CO_2-eq m^{-2} s^{-1})');
title('(c)');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);

subplot(2,3,4);
hold on; box on;
errorbar(bin_v_avg1, v_ch4_avg1a, v_ch4_se1a, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_v_avg2, v_ch4_avg2a, v_ch4_se2a, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_v_avg3, v_ch4_avg3a, v_ch4_se3a, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_v_avg4, v_ch4_avg4a, v_ch4_se4a, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
ylim([0 0.1]);
xlim([0 2.5]);
xlabel('VPD (kPa)');
ylabel('CH_4 flux (mg CO_2-eq m^{-2} s^{-1})');
title('(d)');
set(gca, 'XTick', 0:0.5:2.5, 'FontName', 'Times New Roman', 'FontSize', 12);

subplot(2,3,5);
hold on; box on;
errorbar(bin_Ts_avg1, Ts_ch4_avg1a, Ts_ch4_se1a, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_Ts_avg2, Ts_ch4_avg2a, Ts_ch4_se2a, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_Ts_avg3, Ts_ch4_avg3a, Ts_ch4_se3a, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_Ts_avg4, Ts_ch4_avg4a, Ts_ch4_se4a, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
ylim([-0.05 0.15]);
xlim([-10 30]);
xlabel('T_{soil} 10 cm (\circC)');
ylabel('CH_4 flux (mg CO_2-eq m^{-2} s^{-1})');
title('(e)');
set(gca, 'YTick', -0.05:0.05:0.15, 'FontName', 'Times New Roman', 'FontSize', 12);

subplot(2,3,6);
hold on; box on;
errorbar(bin_sm_avg1, sm_ch4_avg1a, sm_ch4_se1a, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_sm_avg2, sm_ch4_avg2a, sm_ch4_se2a, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_sm_avg3, sm_ch4_avg3a, sm_ch4_se3a, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_sm_avg4, sm_ch4_avg4a, sm_ch4_se4a, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
ylim([0 0.1]);
xlim([0.1 0.6]);
xlabel('SWC 10 cm (m^{3} m^{-3})');
ylabel('CH_4 flux (mg CO_2-eq m^{-2} s^{-1})');
title('(f)');
set(gca, 'XTick', 0.1:0.1:0.6, 'FontName', 'Times New Roman', 'FontSize', 12);


figure(12)
set(gcf, 'Position', [100 100 1200 650]);
subplot(2,3,1);
hold on; box on;
errorbar(bin_Rg_avg1, Rg_ch4_avg1b, Rg_ch4_se1b, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_Rg_avg2, Rg_ch4_avg2b, Rg_ch4_se2b, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_Rg_avg3, Rg_ch4_avg3b, Rg_ch4_se3b, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_Rg_avg4, Rg_ch4_avg4b, Rg_ch4_se4b, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
ylim([0 0.15]);
xlim([0 1500]);
xlabel('R_g (W m^{-2})');
ylabel('CH_4 flux (mg CO_2-eq m^{-2} s^{-1})');
legend('W', 'ST', 'GS', 'SF', 'Location', 'northeast', 'Orientation', 'horizontal', 'FontSize' ,10, 'NumColumns', 2);
title('(a)');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);

subplot(2,3,2);
hold on; box on;
errorbar(bin_Ta_avg1, Ta_ch4_avg1b, Ta_ch4_se1b, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_Ta_avg2, Ta_ch4_avg2b, Ta_ch4_se2b, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_Ta_avg3, Ta_ch4_avg3b, Ta_ch4_se3b, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_Ta_avg4, Ta_ch4_avg4b, Ta_ch4_se4b, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
ylim([-0.05 0.2]);
xlim([-30 30]);
xlabel('T_{air} (\circC)');
ylabel('CH_4 flux (mg CO_2-eq m^{-2} s^{-1})');
title('(b)');
set(gca, 'XTick', -30:10:30, 'FontName', 'Times New Roman', 'FontSize', 12);

subplot(2,3,3);
hold on; box on;
errorbar(bin_p_avg1(~isnan(p_ch4_avg1b)), p_ch4_avg1b(~isnan(p_ch4_avg1b)), p_ch4_se1b(~isnan(p_ch4_avg1b)),...
    'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_p_avg2(~isnan(p_ch4_avg2b)), p_ch4_avg2b(~isnan(p_ch4_avg2b)), p_ch4_se2b(~isnan(p_ch4_avg2b)),...
    'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_p_avg3(~isnan(p_ch4_avg3b)), p_ch4_avg3b(~isnan(p_ch4_avg3b)), p_ch4_se3b(~isnan(p_ch4_avg3b)),...
    'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_p_avg4(~isnan(p_ch4_avg4b)), p_ch4_avg4b(~isnan(p_ch4_avg4b)), p_ch4_se4b(~isnan(p_ch4_avg4b)),...
    'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
ylim([0 0.2]);
xlabel('PPT (mm)');
ylabel('CH_4 flux (mg CO_2-eq m^{-2} s^{-1})');
title('(c)');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);

subplot(2,3,4);
hold on; box on;
errorbar(bin_v_avg1, v_ch4_avg1b, v_ch4_se1b, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_v_avg2, v_ch4_avg2b, v_ch4_se2b, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_v_avg3, v_ch4_avg3b, v_ch4_se3b, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_v_avg4, v_ch4_avg4b, v_ch4_se4b, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
ylim([0 0.2]);
xlim([0 2.5]);
xlabel('VPD (kPa)');
ylabel('CH_4 flux (mg CO_2-eq m^{-2} s^{-1})');
title('(d)');
set(gca, 'XTick', 0:0.5:2.5, 'FontName', 'Times New Roman', 'FontSize', 12);

subplot(2,3,5);
hold on; box on;
errorbar(bin_Ts_avg1, Ts_ch4_avg1b, Ts_ch4_se1b, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_Ts_avg2, Ts_ch4_avg2b, Ts_ch4_se2b, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_Ts_avg3, Ts_ch4_avg3b, Ts_ch4_se3b, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_Ts_avg4, Ts_ch4_avg4b, Ts_ch4_se4b, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
ylim([-0.05 0.25]);
xlim([-10 30]);
xlabel('T_{soil} 10 cm (\circC)');
ylabel('CH_4 flux (mg CO_2-eq m^{-2} s^{-1})');
title('(e)');
set(gca, 'YTick', -0.05:0.05:0.25, 'FontName', 'Times New Roman', 'FontSize', 12);

subplot(2,3,6);
hold on; box on;
errorbar(bin_sm_avg1, sm_ch4_avg1b, sm_ch4_se1b, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_sm_avg2, sm_ch4_avg2b, sm_ch4_se2b, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_sm_avg3, sm_ch4_avg3b, sm_ch4_se3b, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_sm_avg4, sm_ch4_avg4b, sm_ch4_se4b, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
ylim([0 0.2]);
xlim([0.1 0.6]);
xlabel('SWC 10 cm (m^{3} m^{-3})');
ylabel('CH_4 flux (mg CO_2-eq m^{-2} s^{-1})');
title('(f)');
set(gca, 'XTick', 0.1:0.1:0.6, 'FontName', 'Times New Roman', 'FontSize', 12);


x1 = mydat.ch4_flux_filtered_umolm2s(mydat.Rg_Wm2>10 & ~isnan(mydat.ch4_flux_filtered_umolm2s) & ~isnan(mydat.co2_flux_filtered_umolm2s));
y1 = mydat.co2_flux_filtered_umolm2s(mydat.Rg_Wm2>10 & ~isnan(mydat.ch4_flux_filtered_umolm2s) & ~isnan(mydat.co2_flux_filtered_umolm2s));
bn = linspace(-0.1, 0.3, 16);
for i = 1:15
    bin_ch4a{i} = x1(x1 >= bn(i) & x1 < bn(i+1));
    bin_ch4_avga(i,1) = nanmean(bin_ch4a{i});
    bin_co2a{i} = y1(x1 >= bn(i) & x1 < bn(i+1));
    bin_co2_avga(i,1) = nanmean(bin_co2a{i});
    bin_co2_stda(i,1) = nanstd(bin_co2a{i});
    bin_co2_sea(i,1) = bin_co2_stda(i,1)/sqrt(length(bin_co2a{i}));
end

x2 = mydat.ch4_flux_filtered_umolm2s(mydat.Rg_Wm2<=10 & ~isnan(mydat.ch4_flux_filtered_umolm2s) & ~isnan(mydat.co2_flux_filtered_umolm2s));
y2 = mydat.co2_flux_filtered_umolm2s(mydat.Rg_Wm2<=10 & ~isnan(mydat.ch4_flux_filtered_umolm2s) & ~isnan(mydat.co2_flux_filtered_umolm2s));
bn = linspace(-0.2, 0.6, 16);
for i = 1:15
    bin_ch4b{i} = x2(x2 >= bn(i) & x2 < bn(i+1));
    bin_ch4_avgb(i,1) = nanmean(bin_ch4b{i});
    bin_co2b{i} = y2(x2 >= bn(i) & x2 < bn(i+1));
    bin_co2_avgb(i,1) = nanmean(bin_co2b{i});
    bin_co2_stdb(i,1) = nanstd(bin_co2b{i});
    bin_co2_seb(i,1) = bin_co2_stdb(i,1)/sqrt(length(bin_co2b{i}));
end

figure(13);
subplot(1,2,1);
hold on; box on;
plot(x1, y1, '.');
plot(bin_ch4_avga, bin_co2_avga, 'k.', 'MarkerSize', 20);
xlim([-0.2 0.4]);
ylim([-30 20]);
xlabel('CH_4 flux (\mumol m^{-2} s^{-1})');
ylabel('CO_2 flux (\mumol m^{-2} s^{-1})');
title('(a) Daytime');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);

subplot(1,2,2);
hold on; box on;
plot(x2, y2, '.');
plot(bin_ch4_avgb, bin_co2_avgb, 'k.', 'MarkerSize', 20);
xlim([-0.4 0.6]);
ylim([-20 30]);
xlabel('CH_4 flux (\mumol m^{-2} s^{-1})');
ylabel('CO_2 flux (\mumol m^{-2} s^{-1})');
title('(b) Nighttime');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);


figure(14);
hold on; box on;
plot(mydat.TIMESTAMP, cumsum(mydat.co2_flux_filled_umolm2s)*44*1800/1000000);
plot(mydat.TIMESTAMP, cumsum(mydat.NEE_U50_f)*44*1800/1000000);
plot(mydat.TIMESTAMP, cumsum(mydat.NEE_U05_f)*44*1800/1000000);
plot(mydat.TIMESTAMP, cumsum(mydat.NEE_U95_f)*44*1800/1000000);
legend('\itu\rm_* mean', '\itu\rm_* median', '\itu\rm_* low (5%)', '\itu\rm_* high (95%)');
ylabel('Cumulative NEE (g CO_2 m^{-2})');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);


hgs1_ts = timerange('05/01/2014 00:00:00', '08/01/2014 00:00:00');
hgs1 = mydat(hgs1_ts, :);

rg_3 = reshape(hgs1.Rg_Wm2, 48, []);
rg_3_avg = nanmean(rg_3');
rg_3_std = nanstd(rg_3');
rg_3_se = rg_3_std/length(rg_3);

ta_3 = reshape(hgs1.Ta_degC, 48, []);
ta_3_avg = nanmean(ta_3');
ta_3_std = nanstd(ta_3');
ta_3_se = ta_3_std/length(ta_3);

v_3 = reshape(hgs1.VPD_kPa, 48, []);
v_3_avg = nanmean(v_3');
v_3_std = nanstd(v_3');
v_3_se = v_3_std/length(v_3);

ts_3 = reshape(hgs1.Soil_10_degC, 48, []);
ts_3_avg = nanmean(ts_3');
ts_3_std = nanstd(ts_3');
ts_3_se = ts_3_std/length(ts_3);

sm_3 = reshape(hgs1.SWC_10_m3m3,48,[]);
sm_3_avg = nanmean(sm_3');
sm_3_std = nanstd(sm_3');
sm_3_se = sm_3_std/length(sm_3);

co2_3 = reshape(hgs1.co2_flux_filled_umolm2s*44/1000, 48, []);
co2_3_avg = nanmean(co2_3');
co2_3_std = nanstd(co2_3');
co2_3_se = co2_3_std/length(co2_3);

ch4_3 = reshape(hgs1.ch4_flux_filled_umolm2s*16*28/1000, 48, []);
ch4_3_avg = nanmean(ch4_3');
ch4_3_std = nanstd(ch4_3');
ch4_3_se = ch4_3_std/length(ch4_3);

ghg_3 = reshape((hgs1.co2_flux_filled_umolm2s*44/1000+hgs1.ch4_flux_filled_umolm2s*16*28/1000), 48, []);
ghg_3_avg = nanmean(ghg_3');
ghg_3_std = nanstd(ghg_3');
ghg_3_se = ghg_3_std/length(ghg_3);


hgs2_ts = timerange('05/01/2015 00:00:00', '08/01/2015 00:00:00');
hgs2 = mydat(hgs2_ts, :);

rg_4 = reshape(hgs2.Rg_Wm2, 48, []);
rg_4_avg = nanmean(rg_4');
rg_4_std = nanstd(rg_4');
rg_4_se = rg_4_std/length(rg_4);

ta_4 = reshape(hgs2.Ta_degC, 48, []);
ta_4_avg = nanmean(ta_4');
ta_4_std = nanstd(ta_4');
ta_4_se = ta_4_std/length(ta_4);

v_4 = reshape(hgs2.VPD_kPa, 48, []);
v_4_avg = nanmean(v_4');
v_4_std = nanstd(v_4');
v_4_se = v_4_std/length(v_4);

ts_4 = reshape(hgs2.Soil_10_degC, 48, []);
ts_4_avg = nanmean(ts_4');
ts_4_std = nanstd(ts_4');
ts_4_se = ts_4_std/length(ts_4);

sm_4 = reshape(hgs2.SWC_10_m3m3, 48, []);
sm_4_avg = nanmean(sm_4');
sm_4_std = nanstd(sm_4');
sm_4_se = sm_4_std/length(sm_4);

co2_4 = reshape(hgs2.co2_flux_filled_umolm2s*44/1000, 48, []);
co2_4_avg = nanmean(co2_4');
co2_4_std = nanstd(co2_4');
co2_4_se = co2_4_std/length(co2_4);

ch4_4 = reshape(hgs2.ch4_flux_filled_umolm2s*16*28/1000, 48, []);
ch4_4_avg = nanmean(ch4_4');
ch4_4_std = nanstd(ch4_4');
ch4_4_se = ch4_4_std/length(ch4_4);

ghg_4 = reshape((hgs2.co2_flux_filled_umolm2s*44/1000+hgs2.ch4_flux_filled_umolm2s*16*28/1000), 48, []);
ghg_4_avg = nanmean(ghg_4');
ghg_4_std = nanstd(ghg_4');
ghg_4_se = ghg_4_std/length(ghg_4);


hgs3_ts = timerange('05/01/2016 00:00:00', '08/01/2016 00:00:00');
hgs3 = mydat(hgs3_ts, :);

rg_5 = reshape(hgs3.Rg_Wm2, 48 ,[]);
rg_5_avg = nanmean(rg_5');
rg_5_std = nanstd(rg_5');
rg_5_se = rg_5_std/length(rg_5);

ta_5 = reshape(hgs3.Ta_degC, 48, []);
ta_5_avg = nanmean(ta_5');
ta_5_std = nanstd(ta_5');
ta_5_se = ta_5_std/length(ta_5);

v_5 = reshape(hgs3.VPD_kPa, 48, []);
v_5_avg = nanmean(v_5');
v_5_std = nanstd(v_5');
v_5_se = v_5_std/length(v_5);

ts_5 = reshape(hgs3.Soil_10_degC, 48, []);
ts_5_avg = nanmean(ts_5');
ts_5_std = nanstd(ts_5');
ts_5_se = ts_5_std/length(ts_5);

sm_5 = reshape(hgs3.SWC_10_m3m3, 48, []);
sm_5_avg = nanmean(sm_5');
sm_5_std = nanstd(sm_5');
sm_5_se = sm_5_std/length(sm_5);

co2_5 = reshape(hgs3.co2_flux_filled_umolm2s*44/1000, 48, []);
co2_5_avg = nanmean(co2_5');
co2_5_std = nanstd(co2_5');
co2_5_se = co2_5_std/length(co2_5);

ch4_5 = reshape(hgs3.ch4_flux_filled_umolm2s*16*28/1000, 48, []);
ch4_5_avg = nanmean(ch4_5');
ch4_5_std = nanstd(ch4_5');
ch4_5_se = ch4_5_std/length(ch4_5);

ghg_5 = reshape((hgs3.co2_flux_filled_umolm2s*44/1000+hgs3.ch4_flux_filled_umolm2s*16*28/1000), 48, []);
ghg_5_avg = nanmean(ghg_5');
ghg_5_std = nanstd(ghg_5');
ghg_5_se = ghg_5_std/length(ghg_5);


figure(15) %soil thawing
set(gcf, 'Position' ,[20 20 1200 1000]);
subplot(3,3,1);
hold on; box on;
shadedErrorBar(1:48, ghg_3_avg, ghg_3_se, 'lineProps', '-r');
h1 = plot(1:48, ghg_3_avg, 'r', 'LineWidth', 1.5);
shadedErrorBar(1:48, ghg_4_avg, ghg_4_se, 'lineProps', '-b');
h2 = plot(1:48, ghg_4_avg, 'b', 'LineWidth', 1.5);
shadedErrorBar(1:48, ghg_5_avg, ghg_5_se, 'lineProps', '-g');
h3 = plot(1:48, ghg_5_avg, 'g', 'LineWidth', 1.5);
yline(0, '--k');
legend([h1 h2 h3], {'2014','2015','2016'}, 'Location', 'northwest', 'FontSize',9);
xlim([0 48]);
ylim([-0.1 0.2])
xlabel('Hour');
ylabel('net CO_2-eq flux (mg CO_2-eq m^{-2} s^{-1})');
title('(a)');
set(gca, 'XTick', 0:12:48, 'XTickLabel', {'0','6','12','18','24'}, 'FontName', 'Times New Roman', 'FontSize', 14);

subplot(3,3,2);
hold on; box on;
shadedErrorBar(1:48, co2_3_avg, co2_3_se, 'lineProps', '-r');
h1 = plot(1:48, co2_3_avg, 'r', 'LineWidth', 1.5);
shadedErrorBar(1:48, co2_4_avg, co2_4_se, 'lineProps', '-b');
h2 = plot(1:48, co2_4_avg, 'b', 'LineWidth', 1.5);
shadedErrorBar(1:48, co2_5_avg, co2_5_se, 'lineProps', '-g');
h3 = plot(1:48, co2_5_avg, 'g', 'LineWidth', 1.5);
yline(0, '--k');
legend([h1 h2 h3], {'2014','2015','2016'}, 'Location', 'northwest', 'FontSize', 9);
xlim([0 48]);
ylim([-0.1 0.2])
xlabel('Hour');
ylabel('CO_2 flux (mg CO_2 m^{-2} s^{-1})');
title('(b)');
set(gca, 'XTick', 0:12:48, 'XTickLabel', {'0','6','12','18','24'}, 'FontName', 'Times New Roman', 'FontSize', 14);

subplot(3,3,3);
hold on; box on;
shadedErrorBar(1:48, ch4_3_avg, ch4_3_se, 'lineProps', '-r');
h1 = plot(1:48, ch4_3_avg, 'r', 'LineWidth', 1.5);
shadedErrorBar(1:48, ch4_4_avg, ch4_4_se, 'lineProps', '-b');
h2 = plot(1:48, ch4_4_avg, 'b', 'LineWidth', 1.5);
shadedErrorBar(1:48, ch4_5_avg, ch4_5_se, 'lineProps', '-g');
h3 = plot(1:48, ch4_5_avg, 'g', 'LineWidth', 1.5);
yline(0, '--k');
legend([h1 h2 h3], {'2014','2015','2016'}, 'Location', 'northwest', 'FontSize', 9);
xlim([0 48]);
ylim([0 0.02])
xlabel('Hour');
ylabel('CH_4 flux (mg CO_2-eq m^{-2} s^{-1})');
title('(c)');
set(gca,'XTick', 0:12:48, 'XTickLabel', {'0','6','12','18','24'}, 'FontName', 'Times New Roman', 'FontSize', 14);

subplot(3,3,4);
hold on; box on;
shadedErrorBar(1:48, rg_3_avg, rg_3_se, 'lineProps', '-r');
h1 = plot(1:48, rg_3_avg, 'r', 'LineWidth', 1.5);
shadedErrorBar(1:48, rg_4_avg, rg_4_se, 'lineProps', '-b');
h2 = plot(1:48, rg_4_avg, 'b', 'LineWidth', 1.5);
shadedErrorBar(1:48, rg_5_avg, rg_5_se, 'lineProps', '-g');
h3 = plot(1:48, rg_5_avg, 'g', 'LineWidth', 1.5);
legend([h1 h2 h3], {'2014','2015','2016'}, 'Location', 'northwest', 'FontSize', 9);
xlim([0 48]);
ylim([0 1000]);
xlabel('Hour');
ylabel('R_g (W m^{-2})');
title('(d)');
set(gca, 'XTick', 0:12:48, 'XTickLabel', {'0','6','12','18','24'}, 'FontName', 'Times New Roman', 'FontSize', 14);

subplot(3,3,5);
hold on; box on;
shadedErrorBar(1:48, ta_3_avg, ta_3_se, 'lineProps', '-r');
h1 = plot(1:48, ta_3_avg, 'r', 'LineWidth', 1.5);
shadedErrorBar(1:48, ta_4_avg, ta_4_se, 'lineProps', '-b');
h2 = plot(1:48, ta_4_avg, 'b', 'LineWidth', 1.5);
shadedErrorBar(1:48, ta_5_avg, ta_5_se, 'lineProps', '-g');
h3 = plot(1:48, ta_5_avg, 'g', 'LineWidth', 1.5);
legend([h1 h2 h3], {'2014','2015','2016'}, 'Location', 'northwest', 'FontSize', 9);
xlim([0 48]);
ylim([-5 15]);
xlabel('Hour');
ylabel('T_{air} (\circC)');
title('(e)');
set(gca, 'XTick', 0:12:48, 'XTickLabel', {'0','6','12','18','24'}, 'FontName', 'Times New Roman', 'FontSize', 14);

subplot(3,3,6);
hold on; box on;
h = bar(2014:2016, [nansum(st1.PPT_mm), nansum(st2.PPT_mm), nansum(st3.PPT_mm)], 'BarWidth', 0.4);
h.FaceColor = 'flat';
h.CData(1,:) = [1 0 0];
h.CData(2,:) = [0 0 1]';
h.CData(3,:) = [0 1 0];
ylabel('PPT (mm)');
title('(f)');
set(gca, 'XTick', 2014:2016, 'FontName', 'Times New Roman', 'FontSize', 14);

subplot(3,3,7);
hold on; box on;
shadedErrorBar(1:48, v_3_avg, v_3_se, 'lineProps', '-r');
h1 = plot(1:48, v_3_avg, 'r', 'LineWidth', 1.5);
shadedErrorBar(1:48, v_4_avg, v_4_se, 'lineProps', '-b');
h2 = plot(1:48, v_4_avg, 'b', 'LineWidth', 1.5);
shadedErrorBar(1:48, v_5_avg, v_5_se, 'lineProps', '-g');
h3 = plot(1:48, v_5_avg, 'g', 'LineWidth', 1.5);
legend([h1 h2 h3], {'2014','2015','2016'}, 'Location', 'northwest', 'FontSize', 9);
xlim([0 48]);
ylim([0 1]);
xlabel('Hour');
ylabel('VPD (kPa)');
title('(g)');
set(gca, 'XTick', 0:12:48, 'XTickLabel', {'0','6','12','18','24'}, 'FontName', 'Times New Roman', 'FontSize', 14);

subplot(3,3,8);
hold on; box on;
shadedErrorBar(1:48, ts_3_avg, ts_3_se, 'lineProps', '-r');
h1 = plot(1:48, ts_3_avg, 'r', 'LineWidth', 1.5);
shadedErrorBar(1:48, ts_4_avg, ts_4_se, 'lineProps', '-b');
h2 = plot(1:48, ts_4_avg, 'b', 'LineWidth', 1.5);
shadedErrorBar(1:48, ts_5_avg, ts_5_se, 'lineProps', '-g');
h3 = plot(1:48, ts_5_avg, 'g', 'LineWidth', 1.5);
legend([h1 h2 h3], {'2014','2015','2016'}, 'Location', 'northwest', 'FontSize', 9);
xlim([0 48]);
ylim([0 3]);
xlabel('Hour');
ylabel('T_{soil} 10cm (\circC)');
title('(h)');
set(gca, 'XTick', 0:12:48, 'XTickLabel', {'0','6','12','18','24'}, 'FontName', 'Times New Roman', 'FontSize', 14);

subplot(3,3,9);
hold on; box on;
shadedErrorBar(1:48, sm_3_avg, sm_3_se, 'lineProps', '-r');
h1 = plot(1:48, sm_3_avg, 'r', 'LineWidth', 1.5);
shadedErrorBar(1:48, sm_4_avg, sm_4_se, 'lineProps', '-b');
h2 = plot(1:48, sm_4_avg, 'b', 'LineWidth', 1.5);
shadedErrorBar(1:48, sm_5_avg, sm_5_se, 'lineProps', '-g');
h3 = plot(1:48, sm_5_avg, 'g', 'LineWidth', 1.5);
legend([h1 h2 h3], {'2014','2015','2016'}, 'Location', 'east', 'FontSize', 9);
xlim([0 48]);
ylim([0.28 0.38]);
xlabel('Hour');
ylabel('SWC 10cm (m^3 m^{-3})');
title('(i)');
set(gca, 'XTick', 0:12:48, 'XTickLabel', {'0','6','12','18','24'}, 'FontName', 'Times New Roman', 'FontSize', 14);


figure(16) %winter
set(gcf, 'Position' ,[20 20 1200 1000]);
subplot(3,3,1);
hold on; box on;
shadedErrorBar(1:48, ghg_3_avg, ghg_3_se, 'lineProps', '-r');
h1 = plot(1:48, ghg_3_avg, 'r', 'LineWidth', 1.5);
shadedErrorBar(1:48, ghg_4_avg, ghg_4_se, 'lineProps', '-b');
h2 = plot(1:48, ghg_4_avg, 'b', 'LineWidth', 1.5);
yline(0, '--k');
legend([h1 h2], {'2014','2015'}, 'Location', 'southeast', 'FontSize',9);
xlim([0 48]);
ylim([-0.1 0.1])
xlabel('Hour');
ylabel('net CO_2-eq flux (mg CO_2-eq m^{-2} s^{-1})');
title('(a)');
set(gca, 'XTick', 0:12:48, 'XTickLabel', {'0','6','12','18','24'}, 'FontName', 'Times New Roman', 'FontSize', 14);

subplot(3,3,2);
hold on; box on;
shadedErrorBar(1:48, co2_3_avg, co2_3_se, 'lineProps', '-r');
h1 = plot(1:48, co2_3_avg, 'r', 'LineWidth', 1.5);
shadedErrorBar(1:48, co2_4_avg, co2_4_se, 'lineProps', '-b');
h2 = plot(1:48, co2_4_avg, 'b', 'LineWidth', 1.5);
yline(0, '--k');
legend([h1 h2], {'2014','2015'}, 'Location', 'southeast', 'FontSize', 9);
xlim([0 48]);
ylim([-0.1 0.1])
xlabel('Hour');
ylabel('CO_2 flux (mg CO_2 m^{-2} s^{-1})');
title('(b)');
set(gca, 'XTick', 0:12:48, 'XTickLabel', {'0','6','12','18','24'}, 'FontName', 'Times New Roman', 'FontSize', 14);

subplot(3,3,3);
hold on; box on;
shadedErrorBar(1:48, ch4_3_avg, ch4_3_se, 'lineProps', '-r');
h1 = plot(1:48, ch4_3_avg, 'r', 'LineWidth', 1.5);
shadedErrorBar(1:48, ch4_4_avg, ch4_4_se, 'lineProps', '-b');
h2 = plot(1:48, ch4_4_avg, 'b', 'LineWidth', 1.5);
yline(0, '--k');
legend([h1 h2], {'2014','2015'}, 'Location', 'southeast', 'FontSize', 9);
xlim([0 48]);
ylim([0 0.025])
xlabel('Hour');
ylabel('CH_4 flux (mg CO_2-eq m^{-2} s^{-1})');
title('(c)');
set(gca,'XTick', 0:12:48, 'XTickLabel', {'0','6','12','18','24'}, 'FontName', 'Times New Roman', 'FontSize', 14);

subplot(3,3,4);
hold on; box on;
shadedErrorBar(1:48, rg_3_avg, rg_3_se, 'lineProps', '-r');
h1 = plot(1:48, rg_3_avg, 'r', 'LineWidth', 1.5);
shadedErrorBar(1:48, rg_4_avg, rg_4_se, 'lineProps', '-b');
h2 = plot(1:48, rg_4_avg, 'b', 'LineWidth', 1.5);
legend([h1 h2], {'2014','2015'}, 'Location', 'northwest', 'FontSize', 9);
xlim([0 48]);
ylim([0 800]);
xlabel('Hour');
ylabel('R_g (W m^{-2})');
title('(d)');
set(gca, 'XTick', 0:12:48, 'XTickLabel', {'0','6','12','18','24'}, 'FontName', 'Times New Roman', 'FontSize', 14);

subplot(3,3,5);
hold on; box on;
shadedErrorBar(1:48, ta_3_avg, ta_3_se, 'lineProps', '-r');
h1 = plot(1:48, ta_3_avg, 'r', 'LineWidth', 1.5);
shadedErrorBar(1:48, ta_4_avg, ta_4_se, 'lineProps', '-b');
h2 = plot(1:48, ta_4_avg, 'b', 'LineWidth', 1.5);
legend([h1 h2 h3], {'2014','2015'}, 'Location', 'northwest', 'FontSize', 9);
xlim([0 48]);
ylim([-15 5]);
xlabel('Hour');
ylabel('T_{air} (\circC)');
title('(e)');
set(gca, 'XTick', 0:12:48, 'XTickLabel', {'0','6','12','18','24'}, 'FontName', 'Times New Roman', 'FontSize', 14);

subplot(3,3,6);
hold on; box on;
h = bar(2014:2015, [nansum(w1.PPT_mm), nansum(w2.PPT_mm)], 'BarWidth', 0.4);
h.FaceColor = 'flat';
h.CData(1,:) = [1 0 0];
h.CData(2,:) = [0 0 1]';
ylabel('PPT (mm)');
title('(f)');
set(gca, 'XTick', 2014:2015, 'FontName', 'Times New Roman', 'FontSize', 14);

subplot(3,3,7);
hold on; box on;
shadedErrorBar(1:48, v_3_avg, v_3_se, 'lineProps', '-r');
h1 = plot(1:48, v_3_avg, 'r', 'LineWidth', 1.5);
shadedErrorBar(1:48, v_4_avg, v_4_se, 'lineProps', '-b');
h2 = plot(1:48, v_4_avg, 'b', 'LineWidth', 1.5);
legend([h1 h2], {'2014','2015'}, 'Location', 'northwest', 'FontSize', 9);
xlim([0 48]);
ylim([0 0.8]);
xlabel('Hour');
ylabel('VPD (kPa)');
title('(g)');
set(gca, 'XTick', 0:12:48, 'XTickLabel', {'0','6','12','18','24'}, 'FontName', 'Times New Roman', 'FontSize', 14);

subplot(3,3,8);
hold on; box on;
shadedErrorBar(1:48, ts_3_avg, ts_3_se, 'lineProps', '-r');
h1 = plot(1:48, ts_3_avg, 'r', 'LineWidth', 1.5);
shadedErrorBar(1:48, ts_4_avg, ts_4_se, 'lineProps', '-b');
h2 = plot(1:48, ts_4_avg, 'b', 'LineWidth', 1.5);
legend([h1 h2], {'2014','2015'}, 'Location', 'southeast', 'FontSize', 9);
xlim([0 48]);
ylim([-2.5 -0.5]);
xlabel('Hour');
ylabel('T_{soil} 10cm (\circC)');
title('(h)');
set(gca, 'XTick', 0:12:48, 'XTickLabel', {'0','6','12','18','24'}, 'FontName', 'Times New Roman', 'FontSize', 14);

subplot(3,3,9);
hold on; box on;
shadedErrorBar(1:48, sm_3_avg, sm_3_se, 'lineProps', '-r');
h1 = plot(1:48, sm_3_avg, 'r', 'LineWidth', 1.5);
shadedErrorBar(1:48, sm_4_avg, sm_4_se, 'lineProps', '-b');
h2 = plot(1:48, sm_4_avg, 'b', 'LineWidth', 1.5);
legend([h1 h2], {'2014','2015'}, 'Location', 'southeast', 'FontSize', 9);
xlim([0 48]);
ylim([0.2 0.23]);
xlabel('Hour');
ylabel('SWC 10cm (m^3 m^{-3})');
title('(i)');
set(gca, 'XTick', 0:12:48, 'XTickLabel', {'0','6','12','18','24'}, 'FontName', 'Times New Roman', 'FontSize', 14);


figure(17) %soil freezing
set(gcf, 'Position' ,[20 20 1200 1000]);
subplot(3,3,1);
hold on; box on;
shadedErrorBar(1:48, ghg_3_avg, ghg_3_se, 'lineProps', '-r');
h1 = plot(1:48, ghg_3_avg, 'r', 'LineWidth', 1.5);
shadedErrorBar(1:48, ghg_4_avg, ghg_4_se, 'lineProps', '-b');
h2 = plot(1:48, ghg_4_avg, 'b', 'LineWidth', 1.5);
yline(0, '--k');
legend([h1 h2], {'2014','2015'}, 'Location', 'southeast', 'FontSize',9);
xlim([0 48]);
ylim([-0.1 0.2])
xlabel('Hour');
ylabel('net CO_2-eq flux (mg CO_2-eq m^{-2} s^{-1})');
title('(a)');
set(gca, 'XTick', 0:12:48, 'XTickLabel', {'0','6','12','18','24'}, 'FontName', 'Times New Roman', 'FontSize', 14);

subplot(3,3,2);
hold on; box on;
shadedErrorBar(1:48, co2_3_avg, co2_3_se, 'lineProps', '-r');
h1 = plot(1:48, co2_3_avg, 'r', 'LineWidth', 1.5);
shadedErrorBar(1:48, co2_4_avg, co2_4_se, 'lineProps', '-b');
h2 = plot(1:48, co2_4_avg, 'b', 'LineWidth', 1.5);
yline(0, '--k');
legend([h1 h2], {'2014','2015'}, 'Location', 'southeast', 'FontSize', 9);
xlim([0 48]);
ylim([-0.1 0.2])
xlabel('Hour');
ylabel('CO_2 flux (mg CO_2 m^{-2} s^{-1})');
title('(b)');
set(gca, 'XTick', 0:12:48, 'XTickLabel', {'0','6','12','18','24'}, 'FontName', 'Times New Roman', 'FontSize', 14);

subplot(3,3,3);
hold on; box on;
shadedErrorBar(1:48, ch4_3_avg, ch4_3_se, 'lineProps', '-r');
h1 = plot(1:48, ch4_3_avg, 'r', 'LineWidth', 1.5);
shadedErrorBar(1:48, ch4_4_avg, ch4_4_se, 'lineProps', '-b');
h2 = plot(1:48, ch4_4_avg, 'b', 'LineWidth', 1.5);
yline(0, '--k');
legend([h1 h2], {'2014','2015'}, 'Location', 'southeast', 'FontSize', 9);
xlim([0 48]);
ylim([0.01 0.05])
xlabel('Hour');
ylabel('CH_4 flux (mg CO_2-eq m^{-2} s^{-1})');
title('(c)');
set(gca,'XTick', 0:12:48, 'XTickLabel', {'0','6','12','18','24'}, 'FontName', 'Times New Roman', 'FontSize', 14);

subplot(3,3,4);
hold on; box on;
shadedErrorBar(1:48, rg_3_avg, rg_3_se, 'lineProps', '-r');
h1 = plot(1:48, rg_3_avg, 'r', 'LineWidth', 1.5);
shadedErrorBar(1:48, rg_4_avg, rg_4_se, 'lineProps', '-b');
h2 = plot(1:48, rg_4_avg, 'b', 'LineWidth', 1.5);
legend([h1 h2], {'2014','2015'}, 'Location', 'northwest', 'FontSize', 9);
xlim([0 48]);
ylim([0 800]);
xlabel('Hour');
ylabel('R_g (W m^{-2})');
title('(d)');
set(gca, 'XTick', 0:12:48, 'XTickLabel', {'0','6','12','18','24'}, 'FontName', 'Times New Roman', 'FontSize', 14);

subplot(3,3,5);
hold on; box on;
shadedErrorBar(1:48, ta_3_avg, ta_3_se, 'lineProps', '-r');
h1 = plot(1:48, ta_3_avg, 'r', 'LineWidth', 1.5);
shadedErrorBar(1:48, ta_4_avg, ta_4_se, 'lineProps', '-b');
h2 = plot(1:48, ta_4_avg, 'b', 'LineWidth', 1.5);
legend([h1 h2], {'2014','2015'}, 'Location', 'northwest', 'FontSize', 9);
xlim([0 48]);
ylim([-10 10]);
xlabel('Hour');
ylabel('T_{air} (\circC)');
title('(e)');
set(gca, 'XTick', 0:12:48, 'XTickLabel', {'0','6','12','18','24'}, 'FontName', 'Times New Roman', 'FontSize', 14);

subplot(3,3,6);
hold on; box on;
h = bar(2014:2015, [nansum(sf1.PPT_mm), nansum(sf2.PPT_mm)], 'BarWidth', 0.4);
h.FaceColor = 'flat';
h.CData(1,:) = [1 0 0];
h.CData(2,:) = [0 0 1]';
ylabel('PPT (mm)');
title('(f)');
set(gca, 'XTick', 2014:2015, 'FontName', 'Times New Roman', 'FontSize', 14);

subplot(3,3,7);
hold on; box on;
shadedErrorBar(1:48, v_3_avg, v_3_se, 'lineProps', '-r');
h1 = plot(1:48, v_3_avg, 'r', 'LineWidth', 1.5);
shadedErrorBar(1:48, v_4_avg, v_4_se, 'lineProps', '-b');
h2 = plot(1:48, v_4_avg, 'b', 'LineWidth', 1.5);
legend([h1 h2], {'2014','2015'}, 'Location', 'northwest', 'FontSize', 9);
xlim([0 48]);
ylim([0 1]);
xlabel('Hour');
ylabel('VPD (kPa)');
title('(g)');
set(gca, 'XTick', 0:12:48, 'XTickLabel', {'0','6','12','18','24'}, 'FontName', 'Times New Roman', 'FontSize', 14);

subplot(3,3,8);
hold on; box on;
shadedErrorBar(1:48, ts_3_avg, ts_3_se, 'lineProps', '-r');
h1 = plot(1:48, ts_3_avg, 'r', 'LineWidth', 1.5);
shadedErrorBar(1:48, ts_4_avg, ts_4_se, 'lineProps', '-b');
h2 = plot(1:48, ts_4_avg, 'b', 'LineWidth', 1.5);
legend([h1 h2], {'2014','2015'}, 'Location', 'southeast', 'FontSize', 9);
xlim([0 48]);
ylim([3 5]);
xlabel('Hour');
ylabel('T_{soil} 10cm (\circC)');
title('(h)');
set(gca, 'XTick', 0:12:48, 'XTickLabel', {'0','6','12','18','24'}, 'FontName', 'Times New Roman', 'FontSize', 14);

subplot(3,3,9);
hold on; box on;
shadedErrorBar(1:48, sm_3_avg, sm_3_se, 'lineProps', '-r');
h1 = plot(1:48, sm_3_avg, 'r', 'LineWidth', 1.5);
shadedErrorBar(1:48, sm_4_avg, sm_4_se, 'lineProps', '-b');
h2 = plot(1:48, sm_4_avg, 'b', 'LineWidth', 1.5);
legend([h1 h2], {'2014','2015'}, 'Location', 'southeast', 'FontSize', 9);
xlim([0 48]);
ylim([0.4 0.5]);
xlabel('Hour');
ylabel('SWC 10cm (m^3 m^{-3})');
title('(i)');
set(gca, 'XTick', 0:12:48, 'XTickLabel', {'0','6','12','18','24'}, 'FontName', 'Times New Roman', 'FontSize', 14);


figure(18) %first half of growing seasons
set(gcf, 'Position' ,[20 20 1200 1000]);
subplot(3,3,1);
hold on; box on;
shadedErrorBar(1:48, ghg_3_avg, ghg_3_se, 'lineProps', '-r');
h1 = plot(1:48, ghg_3_avg, 'r', 'LineWidth', 1.5);
shadedErrorBar(1:48, ghg_4_avg, ghg_4_se, 'lineProps', '-b');
h2 = plot(1:48, ghg_4_avg, 'b', 'LineWidth', 1.5);
shadedErrorBar(1:48, ghg_5_avg, ghg_5_se, 'lineProps', '-g');
h3 = plot(1:48, ghg_5_avg, 'g', 'LineWidth', 1.5);
yline(0, '--k');
legend([h1 h2 h3], {'2014','2015','2016'}, 'Location', 'southwest', 'FontSize',9);
xlim([0 48]);
ylim([-0.4 0.4])
xlabel('Hour');
ylabel('net CO_2-eq flux (mg CO_2-eq m^{-2} s^{-1})');
title('(a)');
set(gca, 'XTick', 0:12:48, 'XTickLabel', {'0','6','12','18','24'}, 'FontName', 'Times New Roman', 'FontSize', 14);

subplot(3,3,2);
hold on; box on;
shadedErrorBar(1:48, co2_3_avg, co2_3_se, 'lineProps', '-r');
h1 = plot(1:48, co2_3_avg, 'r', 'LineWidth', 1.5);
shadedErrorBar(1:48, co2_4_avg, co2_4_se, 'lineProps', '-b');
h2 = plot(1:48, co2_4_avg, 'b', 'LineWidth', 1.5);
shadedErrorBar(1:48, co2_5_avg, co2_5_se, 'lineProps', '-g');
h3 = plot(1:48, co2_5_avg, 'g', 'LineWidth', 1.5);
yline(0, '--k');
legend([h1 h2 h3], {'2014','2015','2016'}, 'Location', 'northwest', 'FontSize', 9);
xlim([0 48]);
ylim([-0.4 0.4])
xlabel('Hour');
ylabel('CO_2 flux (mg CO_2 m^{-2} s^{-1})');
title('(b)');
set(gca, 'XTick', 0:12:48, 'XTickLabel', {'0','6','12','18','24'}, 'FontName', 'Times New Roman', 'FontSize', 14);

subplot(3,3,3);
hold on; box on;
shadedErrorBar(1:48, ch4_3_avg, ch4_3_se, 'lineProps', '-r');
h1 = plot(1:48, ch4_3_avg, 'r', 'LineWidth', 1.5);
shadedErrorBar(1:48, ch4_4_avg, ch4_4_se, 'lineProps', '-b');
h2 = plot(1:48, ch4_4_avg, 'b', 'LineWidth', 1.5);
shadedErrorBar(1:48, ch4_5_avg, ch4_5_se, 'lineProps', '-g');
h3 = plot(1:48, ch4_5_avg, 'g', 'LineWidth', 1.5);
yline(0, '--k');
legend([h1 h2 h3], {'2014','2015','2016'}, 'Location', 'southeast', 'FontSize', 9);
xlim([0 48]);
ylim([0.04 0.08])
xlabel('Hour');
ylabel('CH_4 flux (mg CO_2-eq m^{-2} s^{-1})');
title('(c)');
set(gca,'XTick', 0:12:48, 'XTickLabel', {'0','6','12','18','24'}, 'FontName', 'Times New Roman', 'FontSize', 14);

subplot(3,3,4);
hold on; box on;
shadedErrorBar(1:48, rg_3_avg, rg_3_se, 'lineProps', '-r');
h1 = plot(1:48, rg_3_avg, 'r', 'LineWidth', 1.5);
shadedErrorBar(1:48, rg_4_avg, rg_4_se, 'lineProps', '-b');
h2 = plot(1:48, rg_4_avg, 'b', 'LineWidth', 1.5);
shadedErrorBar(1:48, rg_5_avg, rg_5_se, 'lineProps', '-g');
h3 = plot(1:48, rg_5_avg, 'g', 'LineWidth', 1.5);
legend([h1 h2 h3], {'2014','2015','2016'}, 'Location', 'northwest', 'FontSize', 9);
xlim([0 48]);
ylim([0 800]);
xlabel('Hour');
ylabel('R_g (W m^{-2})');
title('(d)');
set(gca, 'XTick', 0:12:48, 'XTickLabel', {'0','6','12','18','24'}, 'FontName', 'Times New Roman', 'FontSize', 14);

subplot(3,3,5);
hold on; box on;
shadedErrorBar(1:48, ta_3_avg, ta_3_se, 'lineProps', '-r');
h1 = plot(1:48, ta_3_avg, 'r', 'LineWidth', 1.5);
shadedErrorBar(1:48, ta_4_avg, ta_4_se, 'lineProps', '-b');
h2 = plot(1:48, ta_4_avg, 'b', 'LineWidth', 1.5);
shadedErrorBar(1:48, ta_5_avg, ta_5_se, 'lineProps', '-g');
h3 = plot(1:48, ta_5_avg, 'g', 'LineWidth', 1.5);
legend([h1 h2 h3], {'2014','2015','2016'}, 'Location', 'northwest', 'FontSize', 9);
xlim([0 48]);
ylim([2 16]);
xlabel('Hour');
ylabel('T_{air} (\circC)');
title('(e)');
set(gca, 'YTick', 2:2:16, 'XTick', 0:12:48, 'XTickLabel', {'0','6','12','18','24'}, 'FontName', 'Times New Roman', 'FontSize', 14);

subplot(3,3,6);
hold on; box on;
h = bar(2014:2016, [nansum(hgs1.PPT_mm), nansum(hgs2.PPT_mm), nansum(hgs3.PPT_mm)], 'BarWidth', 0.4);
h.FaceColor = 'flat';
h.CData(1,:) = [1 0 0];
h.CData(2,:) = [0 0 1]';
h.CData(3,:) = [0 1 0];
ylabel('PPT (mm)');
title('(f)');
set(gca, 'XTick', 2014:2016, 'FontName', 'Times New Roman', 'FontSize', 14);

subplot(3,3,7);
hold on; box on;
shadedErrorBar(1:48, v_3_avg, v_3_se, 'lineProps', '-r');
h1 = plot(1:48, v_3_avg, 'r', 'LineWidth', 1.5);
shadedErrorBar(1:48, v_4_avg, v_4_se, 'lineProps', '-b');
h2 = plot(1:48, v_4_avg, 'b', 'LineWidth', 1.5);
shadedErrorBar(1:48, v_5_avg, v_5_se, 'lineProps', '-g');
h3 = plot(1:48, v_5_avg, 'g', 'LineWidth', 1.5);
legend([h1 h2 h3], {'2014','2015','2016'}, 'Location', 'northwest', 'FontSize', 9);
xlim([0 48]);
ylim([0 1]);
xlabel('Hour');
ylabel('VPD (kPa)');
title('(g)');
set(gca, 'XTick', 0:12:48, 'XTickLabel', {'0','6','12','18','24'}, 'FontName', 'Times New Roman', 'FontSize', 14);

subplot(3,3,8);
hold on; box on;
shadedErrorBar(1:48, ts_3_avg, ts_3_se, 'lineProps', '-r');
h1 = plot(1:48, ts_3_avg, 'r', 'LineWidth', 1.5);
shadedErrorBar(1:48, ts_4_avg, ts_4_se, 'lineProps', '-b');
h2 = plot(1:48, ts_4_avg, 'b', 'LineWidth', 1.5);
shadedErrorBar(1:48, ts_5_avg, ts_5_se, 'lineProps', '-g');
h3 = plot(1:48, ts_5_avg, 'g', 'LineWidth', 1.5);
legend([h1 h2 h3], {'2014','2015','2016'}, 'Location', 'northwest', 'FontSize', 9);
xlim([0 48]);
ylim([8 16]);
xlabel('Hour');
ylabel('T_{soil} 10cm (\circC)');
title('(h)');
set(gca, 'XTick', 0:12:48, 'XTickLabel', {'0','6','12','18','24'}, 'FontName', 'Times New Roman', 'FontSize', 14);

subplot(3,3,9);
hold on; box on;
shadedErrorBar(1:48, sm_3_avg, sm_3_se, 'lineProps', '-r');
h1 = plot(1:48, sm_3_avg, 'r', 'LineWidth', 1.5);
shadedErrorBar(1:48, sm_4_avg, sm_4_se, 'lineProps', '-b');
h2 = plot(1:48, sm_4_avg, 'b', 'LineWidth', 1.5);
shadedErrorBar(1:48, sm_5_avg, sm_5_se, 'lineProps', '-g');
h3 = plot(1:48, sm_5_avg, 'g', 'LineWidth', 1.5);
legend([h1 h2 h3], {'2014','2015','2016'}, 'Location', 'east', 'FontSize', 9);
xlim([0 48]);
ylim([0.45 0.49]);
xlabel('Hour');
ylabel('SWC 10cm (m^3 m^{-3})');
title('(i)');
set(gca, 'XTick', 0:12:48, 'XTickLabel', {'0','6','12','18','24'}, 'FontName', 'Times New Roman', 'FontSize', 14);

