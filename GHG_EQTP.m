%Eastern TQP GHG study
%J. Chi
%Jul 15, 2020

clear all; clc; close all;
warning off;
mydat = readtable('C:/Users/jihi0001/Dropbox/Publications/Peng_Chi_Manuscript/Data/GHG_EQTP_lv1.csv');


%30-min non-gapfilled data
mydat.co2 = mydat.co2_flux_filtered_umolm2s * 44 / 1000; %mg CO2-eq m-2 s-1
mydat.ch4_eq1 = mydat.ch4_flux_filtered_umolm2s * 16 * 28 / 1000; %mg CO2-eq m-2 s-1 for GWP28
mydat.ch4_eq2 = mydat.ch4_flux_filtered_umolm2s * 16 * 45 / 1000; % mg CO2-eq m-2 s-1 for SGWP45
mydat.ghg1 = mydat.co2 + mydat.ch4_eq1;
mydat.ghg2 = mydat.co2 + mydat.ch4_eq2;

figure(1)
hold on;
plot(mydat.TIMESTAMP, mydat.co2_flux_filled_umolm2s, '.');
plot(mydat.TIMESTAMP, mydat.co2_flux_filtered_umolm2s, '.');


figure(2)
hold on;
plot(mydat.TIMESTAMP, mydat.ch4_flux_filled_umolm2s, '.');
plot(mydat.TIMESTAMP, mydat.ch4_flux_filtered_umolm2s, '.');

ts = unique(datenum(datestr(mydat.TIMESTAMP, 'yyyy/mm/dd')));
datno = ts(1:end-1);
datstr = datetime(datno, 'ConvertFrom', 'datenum', 'Format', 'yyyy/MM/dd');

%daily mean or sum data
Rg = mean(reshape(mydat.Rg_Wm2, 48,[]));
Ta = mean(reshape(mydat.Ta_degC, 48, []));
VPD = mean(reshape(mydat.VPD_kPa, 48, []));
PPT = sum(reshape(mydat.PPT_mm, 48, []));
Ts10 = mean(reshape(mydat.Soil_10_degC, 48, []));
Ts25 = mean(reshape(mydat.Soil_25_degC, 48, []));
Ts40 = mean(reshape(mydat.Soil_40_degC, 48, []));
SWC = mean(reshape(mydat.SWC_10_m3m3, 48, []));
co2_d = sum(reshape(mydat.co2_flux_filled_umolm2s, 48, [])) * 44 * 1800 / 1000000; %unit in g co2 m-2 d-1
ch4_d1 = sum(reshape(mydat.ch4_flux_filled_umolm2s, 48, [])) * 16 * 28 * 1800 / 1000000; %unit in g co2 m-2 d-1
ch4_d2 = sum(reshape(mydat.ch4_flux_filled_umolm2s, 48, [])) * 16 * 45 * 1800 / 1000000; %unit in g co2 m-2 d-1
ghg_d1 = co2_d + ch4_d1;
ghg_d2 = co2_d + ch4_d2;

%monthly mean or sum data
month1 = datevec(mydat.TIMESTAMP);
monthw(1, :) = [2013, 12, 1, 0, 0, 0];
monthw(2:length(month1), :) = month1(1:end-1, :);
a = unique(month1(:, 1:2), 'rows');
mon = datetime([a, ones([length(a), 1]), zeros([length(a), 1]), zeros([length(a), 1]), zeros([length(a), 1])]);
[b, c, c] = unique(monthw(:, 1:2), 'rows');

co2_monthly = [b, accumarray(c, mydat.co2_flux_filled_umolm2s, [], @sum)* 44 * 1800 / 1000000]; %unit in t co2-eq ha-1 month-1
ch4_monthly1 = [b, accumarray(c, mydat.ch4_flux_filled_umolm2s, [], @sum) * 16 * 28 * 1800 / 1000000]; %unit in t co2-eq ha-1 month-1
ch4_monthly2 = [b, accumarray(c, mydat.ch4_flux_filled_umolm2s, [], @sum) * 16 * 45 * 1800 / 1000000]; 
ghg_monthly1 = co2_monthly(:, end) + ch4_monthly1(:, end);
ghg_monthly2 = co2_monthly(:, end) + ch4_monthly2(:, end);
Ta_monthly = [b, accumarray(c, mydat.Ta_degC, [], @mean)];
Ts10_monthly = [b, accumarray(c, mydat.Soil_10_degC, [], @mean)];
Ts25_monthly = [b, accumarray(c, mydat.Soil_25_degC, [], @mean)];
Ts40_monthly = [b, accumarray(c, mydat.Soil_40_degC, [], @mean)];
VPD_monthly = [b, accumarray(c, mydat.VPD_kPa, [], @mean)];
PPT_monthly = [b, accumarray(c, mydat.PPT_mm, [], @sum)];
Rg_monthly = [b, accumarray(c, mydat.Rg_Wm2, [], @mean)];
SWC_monthly = [b, accumarray(c, mydat.SWC_10_m3m3, [], @mean)];


figure(3) %daily met plots
set(gcf, 'Position', [100 100 800 800]);
subplot(3, 1, 1);
hold on; box on;
yyaxis left
h1 = bar(datstr, Rg, 'FaceColor', [0.8, 0.08, 0.08]);
xline(datetime('2014-01-01 00:30:00'));
xline(datetime('2015-01-01 00:30:00'));
xline(datetime('2016-01-01 00:30:00'));
ylim([0 500]);
ylabel('R_g (W m^{-2})');
set(gca, 'XTick', mon, 'XTickLabel', {'D', 'J', 'F', 'M', 'A', 'M', 'J', 'J',...
    'A', 'S', 'O', 'N', 'D', 'J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D',...
    'J', 'F', 'M', 'A', 'M', 'J', 'J', 'A'},...
    'YTick', [0 100 200 300 400 500], 'YTickLabel', {'0', '100', '200', '300', '400', '500'});
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14, 'Layer', 'top');
yyaxis right
h2 = plot(datstr, Ta, 'LineWidth', 1.5, 'Color', [0.93, 0.69, 0.13]);
ylim([-40 20]);
ylabel('T_{air} (\circC)');
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14, 'Layer', 'top');

subplot(3, 1, 2);
hold on; box on;
plot(datstr, Ts10, 'LineWidth', 1.5, 'Color', 'k');
plot(datstr, Ts25, 'LineWidth', 1.5, 'Color', [0.47, 0.67, 0.19]);
plot(datstr, Ts40, 'LineWidth', 1.5, 'Color', [0.30, 0.75, 0.93]);
xline(datetime('2014-01-01 00:30:00'));
xline(datetime('2015-01-01 00:30:00'));
xline(datetime('2016-01-01 00:30:00'));
xlim([datstr(1) datstr(end)]);
ylim([-10 30]);
set(gca, 'XTick', mon, 'XTickLabel', {'D', 'J', 'F', 'M', 'A', 'M', 'J', 'J',...
    'A', 'S', 'O', 'N', 'D', 'J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D',...
    'J', 'F', 'M', 'A', 'M', 'J', 'J', 'A'},...
    'YTick', [-10 0 10 20 30], 'YTickLabel', {'-10', '0', '10', '20', '30'});
ylabel('T_{soil} (\circC)');
legend('10 cm', '25 cm', '40 cm', 'Location', 'northwest', 'FontSize', 12);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14, 'Layer', 'top');

subplot(3, 1, 3);
hold on; box on;
yyaxis left
plot(datstr, VPD, 'LineWidth', 1.5, 'Color', [0.49, 0.18, 0.56]);
xline(datetime('2014-01-01 00:30:00'));
xline(datetime('2015-01-01 00:30:00'));
xline(datetime('2016-01-01 00:30:00'));
ylim([0 1.5]);
ylabel('VPD (kPa)');
set(gca, 'XTick', mon, 'XTickLabel', {'D', 'J', 'F', 'M', 'A', 'M', 'J', 'J',...
    'A', 'S', 'O', 'N', 'D', 'J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D',...
    'J', 'F', 'M', 'A', 'M', 'J', 'J', 'A'},...
    'YTick', [0 0.5 1 1.5], 'YTickLabel', {'0', '0.5', '1', '1.5'});
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14, 'Layer', 'top');
yyaxis right
bar(datstr, PPT, 'FaceColor', 'b');
ylim([0 40]);
ylabel('PPT (mm)');
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14, 'Layer', 'top', 'YDir', 'reverse');


figure(4) %monthly met plots
set(gcf, 'Position', [100 100 700 800]);
subplot(3, 1, 1);
hold on; box on;
yyaxis left
h1 = bar(Rg_monthly(:, end), 'FaceColor', [1.00, 0.41, 0.16], 'BarWidth', 0.6);
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
h2 = plot(Ta_monthly(:, end), '-o', 'Color', [0.64, 0.08, 0.18], 'LineWidth', 1.5, 'MarkerFaceColor', [0.64, 0.08, 0.18]);
ylim([-40 20]);
ylabel('T_{air} (\circC)');
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14, 'Layer', 'top');

subplot(3, 1, 2);
hold on; box on;
yyaxis right
plot(2:33, VPD_monthly(:, end), '-o', 'LineWidth', 1.5, 'Color', [0.00, 0.00, 1.00], 'MarkerFaceColor', [0.00, 0.00, 1.00]);
xline(2.5, '--k'); xline(13.5, '--k'); xline(26.5, '--k');
xlim([1 34]); ylim([0 0.5]);
ylabel('VPD (kPa)');
set(gca, 'XTick', 2:33, 'XTickLabel', {'D', 'J', 'F', 'M', 'A', 'M', 'J', 'J',...
    'A', 'S', 'O', 'N', 'D', 'J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D',...
    'J', 'F', 'M', 'A', 'M', 'J', 'J', 'A'}, 'YTick', 0:0.1:0.5);
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
bar(2:33, PPT_monthly(:, end), 'FaceColor', [0.30, 0.75, 0.93], 'BarWidth', 0.6);
ylim([0 300]);
ylabel('PPT (mm)');
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
set(gca, 'YTick', 0:50:300);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14, 'Layer', 'top');

subplot(3, 1, 3);
hold on; box on;
yyaxis right
h1 = plot(2:33, Ts10_monthly(:, end), '-o', 'LineWidth', 1.5, 'Color', [0.35,0.50,0.13], 'MarkerFaceColor', [0.35, 0.50, 0.13]);
h2 = plot(2:33, Ts25_monthly(:,end), '-o', 'LineWidth', 1.5, 'Color', [0.50,0.50,0.50], 'MarkerFaceColor', [0.50,0.50,0.50]);
h3 = plot(2:33, Ts40_monthly(:,end), '-o', 'LineWidth', 1.5, 'Color', 'k', 'MarkerFaceColor', 'k');
xline(2.5, '--k'); xline(13.5, '--k'); xline(26.5, '--k');
xlim([1 34]); ylim([-5 20]);
set(gca, 'XTick', 2:33, 'XTickLabel', {'D', 'J', 'F', 'M', 'A', 'M', 'J', 'J',...
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
bar(2:33, SWC_monthly(:, end), 'FaceColor', [0.64, 0.88, 0.46], 'BarWidth', 0.6);
ylim([0 1]);
ylabel('SWC 10 cm (m^{3} m^{-3})');
legend([h1 h2 h3], {'10 cm', '25 cm', '40 cm'}, 'Location', 'northwest', 'FontSize', 12);
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
set(gca, 'YTick', 0:0.2:1);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14, 'Layer', 'top');


figure(5)
set(gcf, 'Position', [100 100 1000 450]);
hold on; box on;
b = bar([co2_monthly(:, end), ch4_monthly1(:, end)], 'BarWidth', 1);
b(1).FaceColor = [0.93, 0.69, 0.13];
plot(ghg_monthly1, 'k-o', 'MarkerFaceColor', 'k');
xline(1.5, '--k');
xline(13.5, '--k');
xline(25.5, '--k');
ylim([-400 300]);
ylabel('CO_2-eq flux (g CO_2-eq m^{-2} month^{-1})');
legend('CO_2', 'CH_4', 'net CO_2-eq', 'Location', 'southeast', 'FontSize', 11);
annotation('textbox', [0.259, 0.86, 0.08, 0.04], 'String', "2014",...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle',...
    'FontName', 'Times New Roman', 'FontSize', 14, 'FontWeight', 'Bold');
annotation('textbox', [0.55, 0.86, 0.08, 0.04], 'String', "2015",...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle',...
    'FontName', 'Times New Roman', 'FontSize', 14, 'FontWeight', 'Bold');
annotation('textbox', [0.78, 0.86, 0.08, 0.04], 'String', "2016",...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle',...
    'FontName', 'Times New Roman', 'FontSize', 14, 'FontWeight', 'Bold');
set(gca, 'YTick', -400:100:300, 'XTick', 1:32, 'XTickLabel', {'D', 'J', 'F', 'M', 'A', 'M', 'J', 'J',...
    'A', 'S', 'O', 'N', 'D', 'J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D',...
    'J', 'F', 'M', 'A', 'M', 'J', 'J', 'A'},...
    'FontName', 'Times New Roman', 'FontSize', 14, 'Layer', 'top');


figure(6)
set(gcf, 'Position', [100 100 1000 450]);
hold on; box on;
b = bar([co2_monthly(:, end), ch4_monthly2(:, end)], 'BarWidth', 1);
b(1).FaceColor = [0.93, 0.69, 0.13];
plot(ghg_monthly2, 'k-o', 'MarkerFaceColor', 'k');
xline(1.5, '--k');
xline(13.5, '--k');
xline(25.5, '--k');
ylim([-400 500]);
ylabel('CO_2-eq flux (g CO_2-eq m^{-2} month^{-1})');
legend('CO_2', 'CH_4', 'net CO_2-eq', 'Location', 'southeast', 'FontSize', 11);
annotation('textbox', [0.259, 0.86, 0.08, 0.04], 'String', "2014",...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle',...
    'FontName', 'Times New Roman', 'FontSize', 14, 'FontWeight', 'Bold');
annotation('textbox', [0.55, 0.86, 0.08, 0.04], 'String', "2015",...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle',...
    'FontName', 'Times New Roman', 'FontSize', 14, 'FontWeight', 'Bold');
annotation('textbox', [0.78, 0.86, 0.08, 0.04], 'String', "2016",...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle',...
    'FontName', 'Times New Roman', 'FontSize', 14, 'FontWeight', 'Bold');
set(gca, 'YTick', -400:100:500, 'XTick', 1:32, 'XTickLabel', {'D', 'J', 'F', 'M', 'A', 'M', 'J', 'J',...
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
pca2 = table(ghg_d1', co2_d', ch4_d1', Ta', Ts10', SWC', Rg', VPD', PPT',...
    'VariableNames', {'GHG', 'CO2', 'CH4', 'Ta', 'Ts10', 'SWC10', 'Rg', 'VPD', 'PPT'});
pca2 = rmmissing(pca2);
pca_norm2 = zscore(table2array(pca2));
[coeff2, ~, latent2, ~, explained2] = pca(pca_norm2);
[R2, p2] = corrcoef(pca_norm2);

%monthly gapfilled data for PCA
pca3 = table(ghg_monthly1, co2_monthly(:, end), ch4_monthly1(:, end),...
    Ta_monthly(:, end), Ts10_monthly(:, end), SWC_monthly(:, end),...
    Rg_monthly(:, end), VPD_monthly(:, end), PPT_monthly(:, end),...
    'VariableNames', {'GHG', 'CO2', 'CH4', 'Ta', 'Ts10', 'SWC10', 'Rg', 'VPD', 'PPT'});
pca3 = rmmissing(pca3);
pca_norm3 = zscore(table2array(pca3));
[coeff3, ~, latent3, ~, explained3] = pca(pca_norm3);
[R3, p3] = corrcoef(pca_norm3);


figure(7)
set(gcf, 'Position', [100 100 1200 350]);
subplot(1, 3, 1);
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
title('Half-hourly');
set(gca, 'XTick', -0.8:0.2:0.8, 'YTick', -0.8:0.2:0.8);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
box on;

subplot(1, 3, 2);
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
title('Daily');
set(gca, 'XTick', -0.8:0.2:0.8, 'YTick', -0.8:0.2:0.8);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
box on;

subplot(1, 3, 3);
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
title('Monthly');
set(gca, 'XTick', -0.8:0.2:0.8, 'YTick', -0.8:0.2:0.8);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
box on;


pca4 = table(mydat.ghg2, mydat.co2, mydat.ch4_eq2,...
    mydat.Ta_degC, mydat.Soil_10_degC, mydat.SWC_10_m3m3, mydat.Rg_Wm2, mydat.VPD_kPa, mydat.PPT_mm, ...
    'VariableNames', {'GHG', 'CO2', 'CH4', 'Ta', 'Ts10', 'SWC10', 'Rg', 'VPD', 'PPT'});
pca4 = rmmissing(pca4);
pca_norm4 = zscore(table2array(pca4));
[coeff4, ~, latent4, ~, explained4] = pca(pca_norm4);
[R4, p4] = corrcoef(pca_norm4);

%daily gapfilled data for PCA
pca5 = table(ghg_d2', co2_d', ch4_d2', Ta', Ts10', SWC', Rg', VPD', PPT',...
    'VariableNames', {'GHG', 'CO2', 'CH4', 'Ta', 'Ts10', 'SWC10', 'Rg', 'VPD', 'PPT'});
pca5 = rmmissing(pca5);
pca_norm5 = zscore(table2array(pca5));
[coeff5, ~, latent5, ~, explained5] = pca(pca_norm5);
[R5, p5] = corrcoef(pca_norm5);

%monthly gapfilled data for PCA
pca6 = table(ghg_monthly2, co2_monthly(:, end), ch4_monthly2(:, end),...
    Ta_monthly(:, end), Ts10_monthly(:, end), SWC_monthly(:, end),...
    Rg_monthly(:, end), VPD_monthly(:, end), PPT_monthly(:, end),...
    'VariableNames', {'GHG', 'CO2', 'CH4', 'Ta', 'Ts10', 'SWC10', 'Rg', 'VPD', 'PPT'});
pca6 = rmmissing(pca6);
pca_norm6 = zscore(table2array(pca6));
[coeff6, ~, latent6, ~, explained6] = pca(pca_norm6);
[R6, p6] = corrcoef(pca_norm6);


figure(8)
set(gcf, 'Position', [100 100 1200 350]);
subplot(1, 3, 1);
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
title('Half-hourly');
set(gca, 'XTick', -0.8:0.2:0.8, 'YTick', -0.8:0.2:0.8);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
box on;

subplot(1, 3, 2);
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
title('Daily');
set(gca, 'XTick', -0.8:0.2:0.8, 'YTick', -0.8:0.2:0.8);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
box on;

subplot(1, 3, 3);
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
title('Monthly');
set(gca, 'XTick', -0.8:0.2:0.8, 'YTick', -0.8:0.2:0.8);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
box on;


[X,Y] = meshgrid(1:9, 1:9);
mycolormap = customcolormap(linspace(0,1,2), [0.93 0.98 1; 0 0.45 0.74]); %blue gradient

% mycolormap = customcolormap(linspace(0,1,3), [0.64 0.08 0.18; 1 1 1; 0 0.45 0.74]); %blue gradient

figure(9)
set(gcf, 'Position', [100 100 1400 820]);
subplot(2, 3, 1);
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
title('Half-hourly');
set(gca, 'XTick', -0.8:0.2:0.8, 'YTick', -0.8:0.2:0.8);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
box on;

subplot(2, 3, 2);
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
title('Daily');
set(gca, 'XTick', -0.8:0.2:0.8, 'YTick', -0.8:0.2:0.8);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
box on;

subplot(2, 3, 3);
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
title('Monthly');
set(gca, 'XTick', -0.8:0.2:0.8, 'YTick', -0.8:0.2:0.8);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
box on;

subplot(2, 3, 4);
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

subplot(2, 3, 5);
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

subplot(2, 3, 6);
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


figure(10)
set(gcf, 'Position', [100 100 1400 820]);
subplot(2, 3, 1);
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
title('Half-hourly');
set(gca, 'XTick', -0.8:0.2:0.8, 'YTick', -0.8:0.2:0.8);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
box on;

subplot(2, 3, 2);
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
title('Daily');
set(gca, 'XTick', -0.8:0.2:0.8, 'YTick', -0.8:0.2:0.8);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
box on;

subplot(2, 3, 3);
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
title('Monthly');
set(gca, 'XTick', -0.8:0.2:0.8, 'YTick', -0.8:0.2:0.8);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
box on;

subplot(2, 3, 4);
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

subplot(2, 3, 5);
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

subplot(2, 3, 6);
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


a1 = find(ismember(mydat.TIMESTAMP, '01/01/2014 00:30')); %update the beginning time
b1 = find(ismember(mydat.TIMESTAMP, '01/01/2015 00:00')); %update the ending time
yr1 = mydat(a1:b1, :); %annual period 2014

c1 = find(ismember(mydat.TIMESTAMP, '03/06/2014 00:30')); %update the beginning time
d1 = find(ismember(mydat.TIMESTAMP, '04/13/2014 00:00')); %update the ending time
st1 = mydat(c1:d1, :); %soil thawing 2014

e1 = find(ismember(mydat.TIMESTAMP, '04/13/2014 00:30')); %update the beginning time
f1 = find(ismember(mydat.TIMESTAMP, '10/02/2014 00:00')); %update the ending time
gs1 = mydat(e1:f1, :); %growing season 2014

g1 = find(ismember(mydat.TIMESTAMP, '10/02/2014 00:30')); %update the beginning time
h1 = find(ismember(mydat.TIMESTAMP, '12/03/2014 00:00')); %update the ending time
sf1  =mydat(g1:h1, :); %soil freezing 2014

w1 = [mydat(a1:c1-1, :); mydat(h1+1:b1, :)];

co2_yr1 = sum(yr1.co2_flux_filled_umolm2s) * 44 * 1800 / 1000000 %unit in g co2-eq m-2 yr-1
ch4_yr1_gwp = sum(yr1.ch4_flux_filled_umolm2s) * 16 * 28 * 1800 / 1000000 %unit in g co2-eq m-2 yr-1
ch4_yr1_sgwp = sum(yr1.ch4_flux_filled_umolm2s) * 16 * 45 * 1800 / 1000000 %unit in g co2-eq m-2 yr-1
ghg_yr1_gwp = co2_yr1 + ch4_yr1_gwp %unit in g co2-eq m-2 yr-1
ghg_yr1_sgwp = co2_yr1 + ch4_yr1_sgwp %unit in g co2-eq m-2 yr-1

co2_gs1 = sum(gs1.co2_flux_filled_umolm2s) * 44 * 1800 / 1000000 %unit in g co2-eq m-2 yr-1
ch4_gs1_gwp = sum(gs1.ch4_flux_filled_umolm2s) * 16 * 28 * 1800/1000000 %unit in g co2-eq m-2 yr-1
ch4_gs1_sgwp = sum(gs1.ch4_flux_filled_umolm2s) * 16 * 45 * 1800/1000000 %unit in g co2-eq m-2 yr-1
ghg_gs1_gwp = co2_gs1 + ch4_gs1_gwp %unit in g co2-eq m-2 yr-1
ghg_gs1_sgwp = co2_gs1 + ch4_gs1_sgwp %unit in g co2-eq m-2 yr-1

co2_sf1 = sum(sf1.co2_flux_filled_umolm2s) * 44 * 1800 / 1000000 %unit in g co2-eq m-2 yr-1
ch4_sf1_gwp = sum(sf1.ch4_flux_filled_umolm2s) * 16 * 28 * 1800 / 1000000 %unit in g co2-eq m-2 yr-1
ch4_sf1_sgwp = sum(sf1.ch4_flux_filled_umolm2s) * 16 * 45 * 1800 / 1000000 %unit in g co2-eq m-2 yr-1
ghg_sf1_gwp = co2_sf1 + ch4_sf1_gwp %unit in g co2-eq m-2 yr-1
ghg_sf1_sgwp = co2_sf1 + ch4_sf1_sgwp %unit in g co2-eq m-2 yr-1

co2_st1 = sum(st1.co2_flux_filled_umolm2s) * 44 * 1800 / 1000000 %unit in g co2-eq m-2 yr-1
ch4_st1_gwp = sum(st1.ch4_flux_filled_umolm2s) * 16 * 28 * 1800 / 1000000 %unit in g co2-eq m-2 yr-1
ch4_st1_sgwp = sum(st1.ch4_flux_filled_umolm2s) * 16 * 45 * 1800 / 1000000 %unit in g co2-eq m-2 yr-1
ghg_st1_gwp = co2_st1 + ch4_st1_gwp %unit in g co2-eq m-2 yr-1
ghg_st1_sgwp = co2_st1 + ch4_st1_sgwp %unit in g co2-eq m-2 yr-1

co2_w1 =co2_yr1 - co2_gs1 - co2_sf1 - co2_st1
ch4_w1_gwp = ch4_yr1_gwp - ch4_gs1_gwp - ch4_sf1_gwp - ch4_st1_gwp
ch4_w1_sgwp = ch4_yr1_sgwp - ch4_gs1_sgwp - ch4_sf1_sgwp - ch4_st1_sgwp
ghg_w1_gwp = co2_w1 + ch4_w1_gwp
ghg_w1_sgwp = co2_w1 + ch4_w1_sgwp


a2 = find(ismember(mydat.TIMESTAMP, '01/01/2015 00:30')); %update the beginning time
b2 = find(ismember(mydat.TIMESTAMP, '01/01/2016 00:00')); %update the ending time
yr2 = mydat(a2:b2, :);

c2 = find(ismember(mydat.TIMESTAMP, '03/14/2015 00:30')); %update the beginning time
d2 = find(ismember(mydat.TIMESTAMP, '05/12/2015 00:00')); %update the ending time
st2 = mydat(c2:d2, :);

e2 = find(ismember(mydat.TIMESTAMP, '05/12/2015 00:30')); %update the beginning time
f2 = find(ismember(mydat.TIMESTAMP, '10/08/2015 00:00')); %update the ending time
gs2 = mydat(e2:f2, :);

g2 = find(ismember(mydat.TIMESTAMP, '10/08/2015 00:30')); %update the beginning time
h2 = find(ismember(mydat.TIMESTAMP, '12/17/2015 00:00')); %update the ending time
sf2 = mydat(g2:h2, :);

w2 = [mydat(a2:c2-1, :); mydat(h2+1:b2, :)];

co2_yr2 = sum(yr2.co2_flux_filled_umolm2s) * 44 * 1800 / 1000000 %unit in g co2-eq m-2 yr-1
ch4_yr2_gwp = sum(yr2.ch4_flux_filled_umolm2s) * 16 * 28 * 1800 / 1000000 %unit in g co2-eq m-2 yr-1
ch4_yr2_sgwp = sum(yr2.ch4_flux_filled_umolm2s) * 16 * 45 * 1800 / 1000000 %unit in g co2-eq m-2 yr-1
ghg_yr2_gwp = co2_yr2 + ch4_yr2_gwp %unit in t CO2-eq ha-1 yr-1
ghg_yr2_sgwp = co2_yr2 + ch4_yr2_sgwp %unit in t CO2-eq ha-1 yr-1

co2_gs2 = sum(gs2.co2_flux_filled_umolm2s) * 44 * 1800 / 1000000 %unit in g co2-eq m-2 yr-1
ch4_gs2_gwp = sum(gs2.ch4_flux_filled_umolm2s) * 16 * 28 * 1800 / 1000000 %unit in g co2-eq m-2 yr-1
ch4_gs2_sgwp = sum(gs2.ch4_flux_filled_umolm2s) * 16 * 45 * 1800 / 1000000 %unit in g co2-eq m-2 yr-1
ghg_gs2_gwp = co2_gs2 + ch4_gs2_gwp %unit in g co2-eq m-2 yr-1
ghg_gs2_sgwp = co2_gs2 + ch4_gs2_sgwp %unit in g co2-eq m-2 yr-1

co2_sf2 = sum(sf2.co2_flux_filled_umolm2s) * 44 * 1800 / 1000000 %unit in g co2-eq m-2 yr-1
ch4_sf2_gwp = sum(sf2.ch4_flux_filled_umolm2s) * 16 * 28 * 1800 / 1000000 %unit in g co2-eq m-2 yr-1
ch4_sf2_sgwp = sum(sf2.ch4_flux_filled_umolm2s) * 16 * 45 * 1800 / 1000000 %unit in g co2-eq m-2 yr-1
ghg_sf2_gwp = co2_sf2 + ch4_sf2_gwp %unit in g co2-eq m-2 yr-1
ghg_sf2_sgwp = co2_sf2 + ch4_sf2_sgwp %unit in g co2-eq m-2 yr-1

co2_st2 = sum(st2.co2_flux_filled_umolm2s) * 44 * 1800 / 1000000 %unit in g co2-eq m-2 yr-1
ch4_st2_gwp = sum(st2.ch4_flux_filled_umolm2s) * 16 * 28 * 1800 / 1000000 %unit in g co2-eq m-2 yr-1
ch4_st2_sgwp = sum(st2.ch4_flux_filled_umolm2s) * 16 * 45 * 1800 / 1000000 %unit in g co2-eq m-2 yr-1
ghg_st2_gwp = co2_st2 + ch4_st2_gwp %unit in g co2-eq m-2 yr-1
ghg_st2_sgwp = co2_st2 + ch4_st2_sgwp %unit in g co2-eq m-2 yr-1

co2_w2 = co2_yr2 - co2_gs2 - co2_sf2 - co2_st2
ch4_w2_gwp = ch4_yr2_gwp - ch4_gs2_gwp - ch4_sf2_gwp - ch4_st2_gwp
ch4_w2_sgwp = ch4_yr2_sgwp - ch4_gs2_sgwp - ch4_sf2_sgwp - ch4_st2_sgwp
ghg_w2_gwp = co2_w2 + ch4_w2_gwp
ghg_w2_sgwp = co2_w2 + ch4_w2_sgwp


a3 = find(ismember(mydat.TIMESTAMP, '01/01/2016 00:30')); %update the beginning time
b3 = find(ismember(mydat.TIMESTAMP, '08/01/2016 00:00')); %update the ending time
yr3 = mydat(a3:b3, :);

c3 = find(ismember(mydat.TIMESTAMP, '03/02/2016 00:30')); %update the beginning time
d3 = find(ismember(mydat.TIMESTAMP, '05/01/2016 00:00')); %update the ending time
st3 = mydat(c3:d3, :);

e3 = find(ismember(mydat.TIMESTAMP, '05/01/2016 00:30')); %update the beginning time
f3 = find(ismember(mydat.TIMESTAMP, '08/01/2016 00:00')); %update the ending time
gs3 = mydat(e3:f3, :);

w3 = mydat(a3:c3-1, :);
w0 = mydat(1:a1-1, :);

w = [w0; w1; w2];
st = [st1; st2; st3];
gs = [gs1; gs2; gs3];
sf = [sf1; sf2];

bn = linspace(0, round(max(w.Rg_Wm2), -2), 16);
for i = 1:15
    bin_Rg1{i} = w.Rg_Wm2(w.Rg_Wm2 >= bn(i) & w.Rg_Wm2 < bn(i+1));
    bin_Rg_avg1(i,1) = nanmean(bin_Rg1{i});
    Rg_co2_1{i} = w.co2(w.Rg_Wm2 >= bn(i) & w.Rg_Wm2 < bn(i+1));
    Rg_ch4_1{i} = w.ch4_eq2(w.Rg_Wm2 >= bn(i) & w.Rg_Wm2 < bn(i+1));
    Rg_ghg1{i} = w.ghg2(w.Rg_Wm2 >= bn(i) & w.Rg_Wm2 < bn(i+1));
    Rg_co2_avg1(i,1) = nanmean(Rg_co2_1{i});
    Rg_ch4_avg1(i,1) = nanmean(Rg_ch4_1{i});
    Rg_ghg_avg1(i,1) = nanmean(Rg_ghg1{i});
    Rg_co2_std1(i,1) = nanstd(Rg_co2_1{i});
    Rg_ch4_std1(i,1) = nanstd(Rg_ch4_1{i});
    Rg_ghg_std1(i,1) = nanstd(Rg_ghg1{i});
    Rg_co2_se1(i,1) = Rg_co2_std1(i,1)/sqrt(length(Rg_co2_1{i}));
    Rg_ch4_se1(i,1) = Rg_ch4_std1(i,1)/sqrt(length(Rg_ch4_1{i}));
    Rg_ghg_se1(i,1) = Rg_ghg_std1(i,1)/sqrt(length(Rg_ghg1{i}));
end

bn = linspace(0, round(max(st.Rg_Wm2), -2), 16);
for i = 1:15
    bin_Rg2{i} = st.Rg_Wm2(st.Rg_Wm2 >= bn(i) & st.Rg_Wm2 < bn(i+1));
    bin_Rg_avg2(i,1) = nanmean(bin_Rg2{i});
    Rg_co2_2{i} = st.co2(st.Rg_Wm2 >= bn(i) & st.Rg_Wm2 < bn(i+1));
    Rg_ch4_2{i} = st.ch4_eq2(st.Rg_Wm2 >= bn(i) & st.Rg_Wm2 < bn(i+1));
    Rg_ghg2{i} = st.ghg2(st.Rg_Wm2 >= bn(i) & st.Rg_Wm2 < bn(i+1));
    Rg_co2_avg2(i,1) = nanmean(Rg_co2_2{i});
    Rg_ch4_avg2(i,1) = nanmean(Rg_ch4_2{i});
    Rg_ghg_avg2(i,1) = nanmean(Rg_ghg2{i});
    Rg_co2_std2(i,1) = nanstd(Rg_co2_2{i});
    Rg_ch4_std2(i,1) = nanstd(Rg_ch4_2{i});
    Rg_ghg_std2(i,1) = nanstd(Rg_ghg2{i});
    Rg_co2_se2(i,1) = Rg_co2_std2(i,1)/sqrt(length(Rg_co2_2{i}));
    Rg_ch4_se2(i,1) = Rg_ch4_std2(i,1)/sqrt(length(Rg_ch4_2{i}));
    Rg_ghg_se2(i,1) = Rg_ghg_std2(i,1)/sqrt(length(Rg_ghg2{i}));
end

bn = linspace(0, round(max(gs.Rg_Wm2), -2), 16);
for i = 1:15
    bin_Rg3{i} = gs.Rg_Wm2(gs.Rg_Wm2 >= bn(i) & gs.Rg_Wm2 < bn(i+1));
    bin_Rg_avg3(i,1) = nanmean(bin_Rg3{i});
    Rg_co2_3{i} = gs.co2(gs.Rg_Wm2 >= bn(i) & gs.Rg_Wm2 < bn(i+1));
    Rg_ch4_3{i} = gs.ch4_eq2(gs.Rg_Wm2 >= bn(i) & gs.Rg_Wm2 < bn(i+1));
    Rg_ghg3{i} = gs.ghg2(gs.Rg_Wm2 >= bn(i) & gs.Rg_Wm2 < bn(i+1));
    Rg_co2_avg3(i,1) = nanmean(Rg_co2_3{i});
    Rg_ch4_avg3(i,1) = nanmean(Rg_ch4_3{i});
    Rg_ghg_avg3(i,1) = nanmean(Rg_ghg3{i});
    Rg_co2_std3(i,1) = nanstd(Rg_co2_3{i});
    Rg_ch4_std3(i,1) = nanstd(Rg_ch4_3{i});
    Rg_ghg_std3(i,1) = nanstd(Rg_ghg3{i});
    Rg_co2_se3(i,1) = Rg_co2_std3(i,1)/sqrt(length(Rg_co2_3{i}));
    Rg_ch4_se3(i,1) = Rg_ch4_std3(i,1)/sqrt(length(Rg_ch4_3{i}));
    Rg_ghg_se3(i,1) = Rg_ghg_std3(i,1)/sqrt(length(Rg_ghg3{i}));
end

bn = linspace(0, round(max(sf.Rg_Wm2), -2), 16);
for i = 1:15
    bin_Rg4{i} = sf.Rg_Wm2(sf.Rg_Wm2 >= bn(i) & sf.Rg_Wm2 < bn(i+1));
    bin_Rg_avg4(i,1) = nanmean(bin_Rg4{i});
    Rg_co2_4{i} = sf.co2(sf.Rg_Wm2 >= bn(i) & sf.Rg_Wm2 < bn(i+1));
    Rg_ch4_4{i} = sf.ch4_eq2(sf.Rg_Wm2 >= bn(i) & sf.Rg_Wm2 < bn(i+1));
    Rg_ghg4{i} = sf.ghg2(sf.Rg_Wm2 >= bn(i) & sf.Rg_Wm2 < bn(i+1));
    Rg_co2_avg4(i,1) = nanmean(Rg_co2_4{i});
    Rg_ch4_avg4(i,1) = nanmean(Rg_ch4_4{i});
    Rg_ghg_avg4(i,1) = nanmean(Rg_ghg4{i});
    Rg_co2_std4(i,1) = nanstd(Rg_co2_4{i});
    Rg_ch4_std4(i,1) = nanstd(Rg_ch4_4{i});
    Rg_ghg_std4(i,1) = nanstd(Rg_ghg4{i});
    Rg_co2_se4(i,1) = Rg_co2_std4(i,1)/sqrt(length(Rg_co2_4{i}));
    Rg_ch4_se4(i,1) = Rg_ch4_std4(i,1)/sqrt(length(Rg_ch4_4{i}));
    Rg_ghg_se4(i,1) = Rg_ghg_std4(i,1)/sqrt(length(Rg_ghg4{i}));
end


figure(11)
set(gcf, 'Position', [100 100 850 800]);
subplot(2, 2, 1);
hold on; box on;
errorbar(bin_Rg_avg1, Rg_ghg_avg1, Rg_ghg_se1,...
    'Marker', 'o', 'MarkerFaceColor', 'k', 'LineWidth', 1.5, 'Color', 'k', 'CapSize', 8, 'MarkerSize', 6);
errorbar(bin_Rg_avg1, Rg_co2_avg1, Rg_co2_se1,...
    'Marker', 'o', 'MarkerFaceColor', [0.93, 0.69, 0.13], 'LineWidth', 1.5, 'Color', [0.93, 0.69, 0.13], 'CapSize', 8, 'MarkerSize', 6);
errorbar(bin_Rg_avg1, Rg_ch4_avg1, Rg_ch4_se1,...
    'Marker', 'o', 'MarkerFaceColor', [0.85, 0.33, 0.10], 'LineWidth', 1.5, 'Color', [0.85, 0.33, 0.10], 'CapSize', 8, 'MarkerSize', 6);
xlim([0 1000]);
ylim([-0.1 0.05]);
xlabel('R_g (W m^{-2})');
ylabel('CO_2-eq flux (mg CO_2-eq m^{-2} s^{-1})');
legend('total CO_2-eq', 'CO_2', 'CH_4', 'Location', 'southwest', 'FontSize', 11);
title('(a) Winter');
set(gca, 'YTick', -0.1:0.05:0.05);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);

subplot(2, 2, 2);
hold on; box on;
errorbar(bin_Rg_avg2, Rg_ghg_avg2, Rg_ghg_se2,...
    'Marker', 'o', 'MarkerFaceColor', 'k', 'LineWidth', 1.5, 'Color', 'k', 'CapSize', 8, 'MarkerSize', 6);
errorbar(bin_Rg_avg2, Rg_co2_avg2, Rg_co2_se2,...
    'Marker', 'o', 'MarkerFaceColor', [0.93, 0.69, 0.13], 'LineWidth', 1.5, 'Color', [0.93, 0.69, 0.13], 'CapSize', 8, 'MarkerSize', 6);
errorbar(bin_Rg_avg2, Rg_ch4_avg2, Rg_ch4_se2,...
    'Marker', 'o', 'MarkerFaceColor', [0.85, 0.33, 0.10], 'LineWidth', 1.5, 'Color', [0.85, 0.33, 0.10], 'CapSize', 8, 'MarkerSize', 6);
xlim([0 1500]);
ylim([-0.1 0.1]);
xlabel('R_g (W m^{-2})');
ylabel('CO_2-eq flux (mg CO_2-eq m^{-2} s^{-1})');
title('(b) Soil thawing');
set(gca, 'YTick', -0.1:0.05:0.1, 'FontName', 'Times New Roman', 'FontSize', 14);

subplot(2, 2, 3);
hold on; box on;
errorbar(bin_Rg_avg3, Rg_ghg_avg3, Rg_ghg_se3,...
    'Marker', 'o', 'MarkerFaceColor', 'k', 'LineWidth', 1.5, 'Color', 'k', 'CapSize', 8, 'MarkerSize', 6);
errorbar(bin_Rg_avg3, Rg_co2_avg3, Rg_co2_se3,...
    'Marker', 'o', 'MarkerFaceColor', [0.93, 0.69, 0.13], 'LineWidth', 1.5, 'Color', [0.93, 0.69, 0.13], 'CapSize', 8, 'MarkerSize', 6);
errorbar(bin_Rg_avg3, Rg_ch4_avg3, Rg_ch4_se3,...
    'Marker', 'o', 'MarkerFaceColor', [0.85, 0.33, 0.10], 'LineWidth', 1.5, 'Color', [0.85, 0.33, 0.10], 'CapSize', 8, 'MarkerSize', 6);
xlim([0 1500]);
ylim([-0.5 0.3]);
xlabel('R_g (W m^{-2})');
ylabel('CO_2-eq flux (mg CO_2-eq m^{-2} s^{-1})');
title('(c) Growing season');
set(gca, 'YTick', -0.5:0.1:0.3, 'FontName', 'Times New Roman', 'FontSize',14);

subplot(2, 2, 4);
hold on; box on;
errorbar(bin_Rg_avg4, Rg_ghg_avg4, Rg_ghg_se4,...
    'Marker', 'o', 'MarkerFaceColor', 'k', 'LineWidth', 1.5, 'Color', 'k', 'CapSize', 8, 'MarkerSize', 6);
errorbar(bin_Rg_avg4, Rg_co2_avg4, Rg_co2_se4,...
    'Marker', 'o', 'MarkerFaceColor', [0.93, 0.69, 0.13], 'LineWidth', 1.5, 'Color', [0.93, 0.69, 0.13], 'CapSize', 8, 'MarkerSize', 6);
errorbar(bin_Rg_avg4, Rg_ch4_avg4, Rg_ch4_se4,...
    'Marker', 'o', 'MarkerFaceColor', [0.85, 0.33, 0.10], 'LineWidth', 1.5, 'Color', [0.85, 0.33, 0.10], 'CapSize', 8, 'MarkerSize', 6);
xlim([0 1500]);
ylim([-0.2 0.15]);
xlabel('R_g (W m^{-2})');
ylabel('CO_2-eq flux (mg CO_2-eq m^{-2} s^{-1})');
title('(d) Soil freezing');
set(gca, 'YTick', -0.2:0.05:0.15, 'FontName', 'Times New Roman', 'FontSize', 14);


bn = linspace(-5.5, 0.5, 16);
for i = 1:15
    bin_Ts1{i} = w.Soil_10_degC(w.Soil_10_degC >= bn(i) & w.Soil_10_degC < bn(i+1));
    bin_Ts_avg1(i,1) = nanmean(bin_Ts1{i});
    Ts_ghg1{i} = w.ghg2(w.Soil_10_degC >= bn(i) & w.Soil_10_degC < bn(i+1));
    Ts_co2_1{i} = w.co2(w.Soil_10_degC >= bn(i) & w.Soil_10_degC < bn(i+1));
    Ts_ch4_1{i} = w.ch4_eq2(w.Soil_10_degC >= bn(i) & w.Soil_10_degC < bn(i+1));
    Ts_ghg_avg1(i,1) = nanmean(Ts_ghg1{i});
    Ts_co2_avg1(i,1) = nanmean(Ts_co2_1{i});
    Ts_ch4_avg1(i,1) = nanmean(Ts_ch4_1{i});
    Ts_ghg_std1(i,1) = nanstd(Ts_ghg1{i});
    Ts_co2_std1(i,1) = nanstd(Ts_co2_1{i});
    Ts_ch4_std1(i,1) = nanstd(Ts_ch4_1{i});
    Ts_ghg_se1(i,1) = Ts_ghg_std1(i,1)/sqrt(length(Ts_ghg1{i}));
    Ts_co2_se1(i,1) = Ts_co2_std1(i,1)/sqrt(length(Ts_co2_1{i}));
    Ts_ch4_se1(i,1) = Ts_ch4_std1(i,1)/sqrt(length(Ts_ch4_1{i}));
end

bn = linspace(-1, 6.5, 16);
for i = 1:15
    bin_Ts2{i} = st.Soil_10_degC(st.Soil_10_degC >= bn(i) & st.Soil_10_degC < bn(i+1));
    bin_Ts_avg2(i,1) = nanmean(bin_Ts2{i});
    Ts_ghg2{i} = st.ghg2(st.Soil_10_degC >= bn(i) & st.Soil_10_degC < bn(i+1));
    Ts_co2_2{i} = st.co2(st.Soil_10_degC >= bn(i) & st.Soil_10_degC < bn(i+1));
    Ts_ch4_2{i} = st.ch4_eq2(st.Soil_10_degC >= bn(i) & st.Soil_10_degC < bn(i+1));
    Ts_ghg_avg2(i,1) = nanmean(Ts_ghg2{i});
    Ts_co2_avg2(i,1) = nanmean(Ts_co2_2{i});
    Ts_ch4_avg2(i,1) = nanmean(Ts_ch4_2{i});
    Ts_ghg_std2(i,1) = nanstd(Ts_ghg2{i});
    Ts_co2_std2(i,1) = nanstd(Ts_co2_2{i});
    Ts_ch4_std2(i,1) = nanstd(Ts_ch4_2{i});
    Ts_ghg_se2(i,1) = Ts_ghg_std2(i,1)/sqrt(length(Ts_ghg2{i}));
    Ts_co2_se2(i,1) = Ts_co2_std2(i,1)/sqrt(length(Ts_co2_2{i}));
    Ts_ch4_se2(i,1)=Ts_ch4_std2(i,1)/sqrt(length(Ts_ch4_2{i}));
end

bn = linspace(0, 23, 16);
for i = 1:15
    bin_Ts3{i} = gs.Soil_10_degC(gs.Soil_10_degC >= bn(i) & gs.Soil_10_degC < bn(i+1));
    bin_Ts_avg3(i,1) = nanmean(bin_Ts3{i});
    Ts_ghg3{i} = gs.ghg2(gs.Soil_10_degC >= bn(i) & gs.Soil_10_degC < bn(i+1));
    Ts_co2_3{i} = gs.co2(gs.Soil_10_degC >= bn(i) & gs.Soil_10_degC < bn(i+1));
    Ts_ch4_3{i} = gs.ch4_eq2(gs.Soil_10_degC >= bn(i) & gs.Soil_10_degC < bn(i+1));
    Ts_ghg_avg3(i,1) = nanmean(Ts_ghg3{i});
    Ts_co2_avg3(i,1) = nanmean(Ts_co2_3{i});
    Ts_ch4_avg3(i,1) = nanmean(Ts_ch4_3{i});
    Ts_ghg_std3(i,1) = nanstd(Ts_ghg3{i});
    Ts_co2_std3(i,1) = nanstd(Ts_co2_3{i});
    Ts_ch4_std3(i,1) = nanstd(Ts_ch4_3{i});
    Ts_ghg_se3(i,1) = Ts_ghg_std3(i,1)/sqrt(length(Ts_ghg3{i}));
    Ts_co2_se3(i,1) = Ts_co2_std3(i,1)/sqrt(length(Ts_co2_3{i}));
    Ts_ch4_se3(i,1) = Ts_ch4_std3(i,1)/sqrt(length(Ts_ch4_3{i}));
end

bn = linspace(-0.5, 12, 16);
for i = 1:15
    bin_Ts4{i} = sf.Soil_10_degC(sf.Soil_10_degC >= bn(i) & sf.Soil_10_degC < bn(i+1));
    bin_Ts_avg4(i,1) = nanmean(bin_Ts4{i});
    Ts_ghg4{i} = sf.ghg2(sf.Soil_10_degC >= bn(i) & sf.Soil_10_degC < bn(i+1));
    Ts_co2_4{i} = sf.co2(sf.Soil_10_degC >= bn(i) & sf.Soil_10_degC < bn(i+1));
    Ts_ch4_4{i} = sf.ch4_eq2(sf.Soil_10_degC >= bn(i) & sf.Soil_10_degC < bn(i+1));
    Ts_ghg_avg4(i,1) = nanmean(Ts_ghg4{i});
    Ts_co2_avg4(i,1) = nanmean(Ts_co2_4{i});
    Ts_ch4_avg4(i,1) = nanmean(Ts_ch4_4{i});
    Ts_ghg_std4(i,1) = nanstd(Ts_ghg4{i});
    Ts_co2_std4(i,1) = nanstd(Ts_co2_4{i});
    Ts_ch4_std4(i,1) = nanstd(Ts_ch4_4{i});
    Ts_ghg_se4(i,1) = Ts_ghg_std4(i,1)/sqrt(length(Ts_ghg4{i}));
    Ts_co2_se4(i,1) = Ts_co2_std4(i,1)/sqrt(length(Ts_co2_4{i}));
    Ts_ch4_se4(i,1) = Ts_ch4_std4(i,1)/sqrt(length(Ts_ch4_4{i}));
end


figure(12)
set(gcf, 'Position', [100 100 850 800]);
subplot(2, 2, 1);
hold on; box on;
errorbar(bin_Ts_avg1, Ts_ghg_avg1, Ts_ghg_se1,...
    'Marker', 'o', 'MarkerFaceColor', 'k', 'LineWidth', 1.5, 'Color', 'k', 'CapSize', 8, 'MarkerSize', 6);
errorbar(bin_Ts_avg1, Ts_co2_avg1, Ts_co2_se1,...
    'Marker', 'o', 'MarkerFaceColor', [0.93, 0.69, 0.13], 'LineWidth', 1.5, 'Color', [0.93, 0.69, 0.13], 'CapSize', 8, 'MarkerSize', 6);
errorbar(bin_Ts_avg1, Ts_ch4_avg1, Ts_ch4_se1,...
    'Marker', 'o', 'MarkerFaceColor', [0.85, 0.33, 0.10], 'LineWidth', 1.5, 'Color', [0.85, 0.33, 0.10], 'CapSize', 8, 'MarkerSize', 6);
xlim([-6 1]);
ylim([-0.04 0.06]);
xlabel('T_{soil} 10 cm (\circC)');
ylabel('CO_2-eq flux (mg CO_2-eq m^{-2} s^{-1})');
legend('net CO_2-eq', 'CO_2', 'CH_4', 'Location', 'southeast', 'FontSize', 11);
title('(a) Winter');
set(gca, 'XTick', -6:1, 'FontName', 'Times New Roman', 'FontSize', 14);

subplot(2, 2, 2);
hold on; box on;
errorbar(bin_Ts_avg2, Ts_ghg_avg2, Ts_ghg_se2,...
    'Marker', 'o', 'MarkerFaceColor', 'k', 'LineWidth', 1.5, 'Color', 'k', 'CapSize', 8, 'MarkerSize', 6);
errorbar(bin_Ts_avg2, Ts_co2_avg2, Ts_co2_se2,...
    'Marker', 'o', 'MarkerFaceColor', [0.93, 0.69, 0.13], 'LineWidth', 1.5, 'Color', [0.93, 0.69, 0.13], 'CapSize', 8, 'MarkerSize', 6);
errorbar(bin_Ts_avg2, Ts_ch4_avg2, Ts_ch4_se2,...
    'Marker', 'o', 'MarkerFaceColor', [0.85, 0.33, 0.10], 'LineWidth', 1.5, 'Color', [0.85, 0.33, 0.10], 'CapSize', 8, 'MarkerSize', 6);
xlim([0 7]);
ylim([-0.02 0.08]);
xlabel('T_{soil} 10 cm (\circC)');
ylabel('CO_2-eq flux (mg CO_2-eq m^{-2} s^{-1})');
title('(b) Soil thawing');
set(gca, 'XTick', 0:7, 'FontName', 'Times New Roman', 'FontSize', 14);

subplot(2, 2, 3);
hold on; box on;
errorbar(bin_Ts_avg3, Ts_ghg_avg3, Ts_ghg_se3,...
    'Marker', 'o', 'MarkerFaceColor', 'k', 'LineWidth', 1.5, 'Color', 'k', 'CapSize', 8, 'MarkerSize', 6);
errorbar(bin_Ts_avg3, Ts_co2_avg3, Ts_co2_se3,...
    'Marker', 'o', 'MarkerFaceColor', [0.93, 0.69, 0.13], 'LineWidth', 1.5, 'Color', [0.93, 0.69, 0.13], 'CapSize', 8, 'MarkerSize', 6);
errorbar(bin_Ts_avg3, Ts_ch4_avg3, Ts_ch4_se3,...
    'Marker', 'o', 'MarkerFaceColor', [0.85, 0.33, 0.10], 'LineWidth', 1.5, 'Color', [0.85, 0.33, 0.10], 'CapSize', 8, 'MarkerSize', 6);
xlim([0 25]);
ylim([-0.3 0.2]);
xlabel('T_{soil} 10 cm (\circC)');
ylabel('CO_2-eq flux (mg CO_2-eq m^{-2} s^{-1})');
title('(c) Growing season');
set(gca,'FontName','Times New Roman','FontSize',14);

subplot(2, 2, 4);
hold on; box on;
errorbar(bin_Ts_avg4, Ts_ghg_avg4, Ts_ghg_se4,...
    'Marker', 'o', 'MarkerFaceColor', 'k', 'LineWidth', 1.5, 'Color', 'k', 'CapSize', 8, 'MarkerSize', 6);
errorbar(bin_Ts_avg4, Ts_co2_avg4, Ts_co2_se4,...
    'Marker', 'o', 'MarkerFaceColor', [0.93, 0.69, 0.13], 'LineWidth', 1.5, 'Color', [0.93, 0.69, 0.13], 'CapSize', 8, 'MarkerSize', 6);
errorbar(bin_Ts_avg4, Ts_ch4_avg4, Ts_ch4_se4,...
    'Marker', 'o', 'MarkerFaceColor', [0.85, 0.33, 0.10], 'LineWidth', 1.5, 'Color', [0.85, 0.33, 0.10], 'CapSize', 8, 'MarkerSize', 6);
xlim([0 12]);
ylim([-0.1 0.15]);
xlabel('T_{soil} 10 cm (\circC)');
ylabel('CO_2-eq flux (mg CO_2-eq m^{-2} s^{-1})');
title('(d) Soil freezing');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);


bn = linspace(-30, 15, 16);
for i = 1:15
    bin_Ta1{i} = w.Ta_degC(w.Ta_degC >= bn(i) & w.Ta_degC < bn(i+1));
    bin_Ta_avg1(i,1) = nanmean(bin_Ta1{i});
    Ta_ghg1{i} = w.ghg2(w.Ta_degC >= bn(i) & w.Ta_degC < bn(i+1));
    Ta_co2_1{i} = w.co2(w.Ta_degC >= bn(i) & w.Ta_degC < bn(i+1));
    Ta_ch4_1{i} = w.ch4_eq2(w.Ta_degC >= bn(i) & w.Ta_degC < bn(i+1));
    Ta_ghg_avg1(i,1) = nanmean(Ta_ghg1{i});
    Ta_co2_avg1(i,1) = nanmean(Ta_co2_1{i}); 
    Ta_ch4_avg1(i,1) = nanmean(Ta_ch4_1{i});
    Ta_ghg_std1(i,1) = nanstd(Ta_ghg1{i});
    Ta_co2_std1(i,1) = nanstd(Ta_co2_1{i});
    Ta_ch4_std1(i,1) = nanstd(Ta_ch4_1{i});
    Ta_ghg_se1(i,1) = Ta_ghg_std1(i,1)/sqrt(length(Ta_ghg1{i}));
    Ta_co2_se1(i,1) = Ta_co2_std1(i,1)/sqrt(length(Ta_co2_1{i}));
    Ta_ch4_se1(i,1) = Ta_ch4_std1(i,1)/sqrt(length(Ta_ch4_1{i}));
end

bn = linspace(-15, 22, 16);
for i = 1:15
    bin_Ta2{i} = st.Ta_degC(st.Ta_degC >= bn(i) & st.Ta_degC < bn(i+1));
    bin_Ta_avg2(i,1) = nanmean(bin_Ta2{i});
    Ta_ghg2{i} = st.ghg2(st.Ta_degC >= bn(i) & st.Ta_degC < bn(i+1));
    Ta_co2_2{i} = st.co2(st.Ta_degC >= bn(i) & st.Ta_degC < bn(i+1));
    Ta_ch4_2{i} = st.ch4_eq2(st.Ta_degC >= bn(i) & st.Ta_degC < bn(i+1));
    Ta_ghg_avg2(i,1) = nanmean(Ta_ghg2{i});
    Ta_co2_avg2(i,1) = nanmean(Ta_co2_2{i});
    Ta_ch4_avg2(i,1) = nanmean(Ta_ch4_2{i});
    Ta_ghg_std2(i,1) = nanstd(Ta_ghg2{i});
    Ta_co2_std2(i,1) = nanstd(Ta_co2_2{i});
    Ta_ch4_std2(i,1) = nanstd(Ta_ch4_2{i});
    Ta_ghg_se2(i,1) = Ta_ghg_std2(i,1)/sqrt(length(Ta_ghg2{i}));
    Ta_co2_se2(i,1) = Ta_co2_std2(i,1)/sqrt(length(Ta_co2_2{i}));
    Ta_ch4_se2(i,1) = Ta_ch4_std2(i,1)/sqrt(length(Ta_ch4_2{i}));
end

bn = linspace(-8, 25, 16);
for i = 1:15
    bin_Ta3{i} = gs.Ta_degC(gs.Ta_degC >= bn(i) & gs.Ta_degC < bn(i+1));
    bin_Ta_avg3(i,1) = nanmean(bin_Ta3{i});
    Ta_ghg3{i} = gs.ghg2(gs.Ta_degC >= bn(i) & gs.Ta_degC < bn(i+1));
    Ta_co2_3{i} = gs.co2(gs.Ta_degC >= bn(i) & gs.Ta_degC < bn(i+1));
    Ta_ch4_3{i} = gs.ch4_eq2(gs.Ta_degC >= bn(i) & gs.Ta_degC < bn(i+1));
    Ta_ghg_avg3(i,1) = nanmean(Ta_ghg3{i});
    Ta_co2_avg3(i,1) = nanmean(Ta_co2_3{i});
    Ta_ch4_avg3(i,1) = nanmean(Ta_ch4_3{i});
    Ta_ghg_std3(i,1) = nanstd(Ta_ghg3{i});
    Ta_co2_std3(i,1) = nanstd(Ta_co2_3{i});
    Ta_ch4_std3(i,1) = nanstd(Ta_ch4_3{i});
    Ta_ghg_se3(i,1) = Ta_ghg_std3(i,1)/sqrt(length(Ta_ghg3{i}));
    Ta_co2_se3(i,1) = Ta_co2_std3(i,1)/sqrt(length(Ta_co2_3{i}));
    Ta_ch4_se3(i,1) = Ta_ch4_std3(i,1)/sqrt(length(Ta_ch4_3{i}));
end

bn = linspace(-20, 20, 16);
for i = 1:15
    bin_Ta4{i} = sf.Ta_degC(sf.Ta_degC >= bn(i) & sf.Ta_degC < bn(i+1));
    bin_Ta_avg4(i,1) = nanmean(bin_Ta4{i});
    Ta_ghg4{i} = sf.ghg2(sf.Ta_degC >= bn(i) & sf.Ta_degC < bn(i+1));
    Ta_co2_4{i} = sf.co2(sf.Ta_degC >= bn(i) & sf.Ta_degC < bn(i+1));
    Ta_ch4_4{i} = sf.ch4_eq2(sf.Ta_degC >= bn(i) & sf.Ta_degC < bn(i+1));
    Ta_ghg_avg4(i,1) = nanmean(Ta_ghg4{i});
    Ta_co2_avg4(i,1) = nanmean(Ta_co2_4{i});
    Ta_ch4_avg4(i,1) = nanmean(Ta_ch4_4{i});
    Ta_ghg_std4(i,1) = nanstd(Ta_ghg4{i});
    Ta_co2_std4(i,1) = nanstd(Ta_co2_4{i});
    Ta_ch4_std4(i,1) = nanstd(Ta_ch4_4{i});
    Ta_ghg_se4(i,1) = Ta_ghg_std4(i,1)/sqrt(length(Ta_ghg4{i}));
    Ta_co2_se4(i,1) = Ta_co2_std4(i,1)/sqrt(length(Ta_co2_4{i}));
    Ta_ch4_se4(i,1) = Ta_ch4_std4(i,1)/sqrt(length(Ta_ch4_4{i}));
end


figure(13)
set(gcf, 'Position', [100 100 850 800]);
subplot(2, 2, 1);
hold on; box on;
errorbar(bin_Ta_avg1, Ta_ghg_avg1, Ta_ghg_se1,...
    'Marker', 'o', 'MarkerFaceColor', 'k', 'LineWidth', 1.5, 'Color', 'k', 'CapSize', 8, 'MarkerSize', 6);
errorbar(bin_Ta_avg1, Ta_co2_avg1, Ta_co2_se1,...
    'Marker', 'o', 'MarkerFaceColor', [0.93, 0.69, 0.13], 'LineWidth', 1.5, 'Color', [0.93, 0.69, 0.13], 'CapSize', 8, 'MarkerSize', 6);
errorbar(bin_Ta_avg1, Ta_ch4_avg1, Ta_ch4_se1,...
    'Marker', 'o', 'MarkerFaceColor', [0.85, 0.33, 0.10], 'LineWidth', 1.5, 'Color', [0.85, 0.33, 0.10], 'CapSize', 8, 'MarkerSize', 6);
xlim([-30 20]);
ylim([-0.15 0.1]);
xlabel('T_{air} (\circC)');
ylabel('CO_2-eq flux (mg CO_2-eq m^{-2} s^{-1})');
legend('net CO_2-eq', 'CO_2', 'CH_4', 'Location', 'southeast', 'FontSize', 11);
title('(a) Winter');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);

subplot(2, 2, 2);
hold on; box on;
errorbar(bin_Ta_avg2, Ta_ghg_avg2, Ta_ghg_se2,...
    'Marker', 'o', 'MarkerFaceColor', 'k', 'LineWidth', 1.5, 'Color', 'k', 'CapSize', 8, 'MarkerSize', 6);
errorbar(bin_Ta_avg2, Ta_co2_avg2, Ta_co2_se2,...
    'Marker', 'o', 'MarkerFaceColor', [0.93, 0.69, 0.13], 'LineWidth', 1.5, 'Color', [0.93, 0.69, 0.13], 'CapSize', 8, 'MarkerSize', 6);
errorbar(bin_Ta_avg2, Ta_ch4_avg2, Ta_ch4_se2,...
    'Marker', 'o', 'MarkerFaceColor', [0.85, 0.33, 0.10], 'LineWidth', 1.5, 'Color', [0.85, 0.33, 0.10], 'CapSize', 8, 'MarkerSize', 6);
xlim([-20 25]);
ylim([-0.04 0.08]);
xlabel('T_{air} (\circC)');
ylabel('CO_2-eq flux (mg CO_2-eq m^{-2} s^{-1})');
title('(b) Soil thawing');
set(gca, 'XTick', -20:5:25, 'FontName', 'Times New Roman', 'FontSize', 14);

subplot(2, 2, 3);
hold on; box on;
errorbar(bin_Ta_avg3, Ta_ghg_avg3, Ta_ghg_se3,...
    'Marker', 'o', 'MarkerFaceColor', 'k', 'LineWidth', 1.5, 'Color', 'k', 'CapSize', 8, 'MarkerSize', 6);
errorbar(bin_Ta_avg3, Ta_co2_avg3, Ta_co2_se3,...
    'Marker', 'o', 'MarkerFaceColor', [0.93, 0.69, 0.13], 'LineWidth', 1.5, 'Color', [0.93, 0.69, 0.13], 'CapSize', 8, 'MarkerSize', 6);
errorbar(bin_Ta_avg3, Ta_ch4_avg3, Ta_ch4_se3,...
    'Marker', 'o', 'MarkerFaceColor', [0.85, 0.33, 0.10], 'LineWidth', 1.5, 'Color', [0.85, 0.33, 0.10], 'CapSize', 8, 'MarkerSize', 6);
xlim([-10 30]);
ylim([-0.4 0.3]);
xlabel('T_{air} (\circC)');
ylabel('CO_2-eq flux (mg CO_2-eq m^{-2} s^{-1})');
title('(c) Growing season');
set(gca, 'YTick', -0.4:0.1:0.3, 'FontName', 'Times New Roman', 'FontSize', 14);

subplot(2, 2, 4);
hold on; box on;
errorbar(bin_Ta_avg4, Ta_ghg_avg4, Ta_ghg_se4,...
    'Marker', 'o', 'MarkerFaceColor', 'k', 'LineWidth', 1.5, 'Color', 'k', 'CapSize', 8, 'MarkerSize', 6);
errorbar(bin_Ta_avg4, Ta_co2_avg4, Ta_co2_se4,...
    'Marker', 'o', 'MarkerFaceColor', [0.93, 0.69, 0.13], 'LineWidth', 1.5, 'Color', [0.93, 0.69, 0.13], 'CapSize', 8, 'MarkerSize', 6);
errorbar(bin_Ta_avg4, Ta_ch4_avg4, Ta_ch4_se4,...
    'Marker', 'o', 'MarkerFaceColor', [0.85, 0.33, 0.10], 'LineWidth', 1.5, 'Color', [0.85, 0.33, 0.10], 'CapSize', 8, 'MarkerSize', 6);
xlim([-25 20]);
ylim([-0.3 0.4]);
xlabel('T_{air} (\circC)');
ylabel('CO_2-eq flux (mg CO_2-eq m^{-2} s^{-1})');
title('(d) Soil freezing');
set(gca, 'YTick', -0.3:0.1:0.4, 'XTick', -25:5:20, 'FontName', 'Times New Roman', 'FontSize', 14);


bn = linspace(0, 1.5, 16);
for i = 1:15
    bin_v1{i} = w.VPD_kPa(w.VPD_kPa >= bn(i) & w.VPD_kPa < bn(i+1));
    bin_v_avg1(i,1) = nanmean(bin_v1{i});
    v_ghg1{i} = w.ghg2(w.VPD_kPa >= bn(i) & w.VPD_kPa < bn(i+1));
    v_co2_1{i} = w.co2(w.VPD_kPa >= bn(i) & w.VPD_kPa < bn(i+1));
    v_ch4_1{i} = w.ch4_eq2(w.VPD_kPa >= bn(i) & w.VPD_kPa < bn(i+1));
    v_ghg_avg1(i,1) = nanmean(v_ghg1{i});
    v_co2_avg1(i,1) = nanmean(v_co2_1{i}); 
    v_ch4_avg1(i,1) = nanmean(v_ch4_1{i});
    v_ghg_std1(i,1) = nanstd(v_ghg1{i});
    v_co2_std1(i,1) = nanstd(v_co2_1{i});
    v_ch4_std1(i,1) = nanstd(v_ch4_1{i});
    v_ghg_se1(i,1) = v_ghg_std1(i,1)/sqrt(length(v_ghg1{i}));
    v_co2_se1(i,1) = v_co2_std1(i,1)/sqrt(length(v_co2_1{i}));
    v_ch4_se1(i,1) = v_ch4_std1(i,1)/sqrt(length(v_ch4_1{i}));
end

bn = linspace(0, 2.3, 16);
for i = 1:15
    bin_v2{i} = st.VPD_kPa(st.VPD_kPa >= bn(i) & st.VPD_kPa < bn(i+1));
    bin_v_avg2(i,1) = nanmean(bin_v2{i});
    v_ghg2{i} = st.ghg2(st.VPD_kPa >= bn(i) & st.VPD_kPa < bn(i+1));
    v_co2_2{i} = st.co2(st.VPD_kPa >= bn(i) & st.VPD_kPa < bn(i+1));
    v_ch4_2{i} = st.ch4_eq2(st.VPD_kPa >= bn(i) & st.VPD_kPa < bn(i+1));
    v_ghg_avg2(i,1) = nanmean(v_ghg2{i});
    v_co2_avg2(i,1) = nanmean(v_co2_2{i});
    v_ch4_avg2(i,1) = nanmean(v_ch4_2{i});
    v_ghg_std2(i,1) = nanstd(v_ghg2{i});
    v_co2_std2(i,1) = nanstd(v_co2_2{i});
    v_ch4_std2(i,1) = nanstd(v_ch4_2{i});
    v_ghg_se2(i,1) = v_ghg_std2(i,1)/sqrt(length(v_ghg2{i}));
    v_co2_se2(i,1) = v_co2_std2(i,1)/sqrt(length(v_co2_2{i}));
    v_ch4_se2(i,1) = v_ch4_std2(i,1)/sqrt(length(v_ch4_2{i}));
end

bn = linspace(0, 2.2, 16);
for i = 1:15
    bin_v3{i} = gs.VPD_kPa(gs.VPD_kPa >= bn(i) & gs.VPD_kPa < bn(i+1));
    bin_v_avg3(i,1) = nanmean(bin_v3{i});
    v_ghg3{i} = gs.ghg2(gs.VPD_kPa >= bn(i) & gs.VPD_kPa < bn(i+1));
    v_co2_3{i} = gs.co2(gs.VPD_kPa >= bn(i) & gs.VPD_kPa < bn(i+1));
    v_ch4_3{i} = gs.ch4_eq2(gs.VPD_kPa >= bn(i) & gs.VPD_kPa < bn(i+1));
    v_ghg_avg3(i,1) = nanmean(v_ghg3{i});
    v_co2_avg3(i,1) = nanmean(v_co2_3{i});
    v_ch4_avg3(i,1) = nanmean(v_ch4_3{i});
    v_ghg_std3(i,1) = nanstd(v_ghg3{i});
    v_co2_std3(i,1) = nanstd(v_co2_3{i});
    v_ch4_std3(i,1) = nanstd(v_ch4_3{i});
    v_ghg_se3(i,1) = v_ghg_std3(i,1)/sqrt(length(v_ghg3{i}));
    v_co2_se3(i,1) = v_co2_std3(i,1)/sqrt(length(v_co2_3{i}));
    v_ch4_se3(i,1) = v_ch4_std3(i,1)/sqrt(length(v_ch4_3{i}));
end

bn = linspace(0, 1.8, 16);
for i = 1:15
    bin_v4{i} = sf.VPD_kPa(sf.VPD_kPa >= bn(i) & sf.VPD_kPa < bn(i+1));
    bin_v_avg4(i,1) = nanmean(bin_v4{i});
    v_ghg4{i} = sf.ghg2(sf.VPD_kPa >= bn(i) & sf.VPD_kPa < bn(i+1));
    v_co2_4{i} = sf.co2(sf.VPD_kPa >= bn(i) & sf.VPD_kPa < bn(i+1));
    v_ch4_4{i} = sf.ch4_eq2(sf.VPD_kPa >= bn(i) & sf.VPD_kPa < bn(i+1));
    v_ghg_avg4(i,1) = nanmean(v_ghg4{i});
    v_co2_avg4(i,1) = nanmean(v_co2_4{i});
    v_ch4_avg4(i,1) = nanmean(v_ch4_4{i});
    v_ghg_std4(i,1) = nanstd(v_ghg4{i});
    v_co2_std4(i,1) = nanstd(v_co2_4{i});
    v_ch4_std4(i,1) = nanstd(v_ch4_4{i});
    v_ghg_se4(i,1) = v_ghg_std4(i,1)/sqrt(length(v_ghg4{i}));
    v_co2_se4(i,1) = v_co2_std4(i,1)/sqrt(length(v_co2_4{i}));
    v_ch4_se4(i,1) = v_ch4_std4(i,1)/sqrt(length(v_ch4_4{i}));
end


figure(14)
set(gcf, 'Position', [100 100 850 800]);
subplot(2, 2, 1);
hold on; box on;
errorbar(bin_v_avg1, v_ghg_avg1, v_ghg_se1,...
    'Marker', 'o', 'MarkerFaceColor', 'k', 'LineWidth', 1.5, 'Color', 'k', 'CapSize', 8, 'MarkerSize', 6);
errorbar(bin_v_avg1, v_co2_avg1, v_co2_se1,...
    'Marker', 'o', 'MarkerFaceColor', [0.93, 0.69, 0.13], 'LineWidth', 1.5, 'Color', [0.93, 0.69, 0.13], 'CapSize', 8, 'MarkerSize', 6);
errorbar(bin_v_avg1, v_ch4_avg1, v_ch4_se1, 'Marker', 'o', 'MarkerFaceColor', [0.85, 0.33, 0.10], 'LineWidth', 1.5, 'Color', [0.85, 0.33, 0.10], 'CapSize', 8, 'MarkerSize', 6);
xlim([0 1.5]);
ylim([-0.06 0.04]);
xlabel('VPD (kPa)');
ylabel('CO_2-eq flux (mg CO_2-eq m^{-2} s^{-1})');
legend('net CO_2-eq', 'CO_2', 'CH_4', 'Location', 'southwest', 'FontSize', 11);
title('(a) Winter');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);

subplot(2, 2, 2);
hold on; box on;
errorbar(bin_v_avg2, v_ghg_avg2, v_ghg_se2,...
    'Marker', 'o', 'MarkerFaceColor', 'k', 'LineWidth', 1.5, 'Color', 'k', 'CapSize', 8, 'MarkerSize', 6);
errorbar(bin_v_avg2, v_co2_avg2, v_co2_se2,...
    'Marker', 'o', 'MarkerFaceColor', [0.93, 0.69, 0.13], 'LineWidth', 1.5, 'Color', [0.93, 0.69, 0.13], 'CapSize', 8, 'MarkerSize', 6);
errorbar(bin_v_avg2, v_ch4_avg2, v_ch4_se2,...
    'Marker', 'o', 'MarkerFaceColor', [0.85, 0.33, 0.10], 'LineWidth', 1.5, 'Color', [0.85, 0.33, 0.10], 'CapSize', 8, 'MarkerSize', 6);
xlim([0 2.5]);
ylim([-0.02 0.1]);
xlabel('VPD (kPa)');
ylabel('CO_2-eq flux (mg CO_2-eq m^{-2} s^{-1})');
title('(b) Soil thawing');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);

subplot(2, 2, 3);
hold on; box on;
errorbar(bin_v_avg3, v_ghg_avg3, v_ghg_se3,...
    'Marker', 'o', 'MarkerFaceColor', 'k', 'LineWidth', 1.5, 'Color', 'k', 'CapSize', 8, 'MarkerSize', 6);
errorbar(bin_v_avg3, v_co2_avg3, v_co2_se3,...
    'Marker', 'o', 'MarkerFaceColor', [0.93, 0.69, 0.13], 'LineWidth', 1.5, 'Color', [0.93, 0.69, 0.13], 'CapSize', 8, 'MarkerSize', 6);
errorbar(bin_v_avg3, v_ch4_avg3, v_ch4_se3,...
    'Marker', 'o', 'MarkerFaceColor', [0.85, 0.33, 0.10], 'LineWidth', 1.5, 'Color',[0.85,0.33,0.10],'CapSize',8,'MarkerSize',6);
xlim([0 2.5]);
ylim([-0.4 0.2]);
xlabel('VPD (kPa)');
ylabel('CO_2-eq flux (mg CO_2-eq m^{-2} s^{-1})');
title('(c) Growing season');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);

subplot(2, 2, 4);
hold on; box on;
errorbar(bin_v_avg4 ,v_ghg_avg4, v_ghg_se4,...
    'Marker', 'o', 'MarkerFaceColor', 'k', 'LineWidth', 1.5, 'Color', 'k', 'CapSize', 8, 'MarkerSize', 6);
errorbar(bin_v_avg4, v_co2_avg4, v_co2_se4,...
    'Marker', 'o', 'MarkerFaceColor', [0.93, 0.69, 0.13], 'LineWidth', 1.5, 'Color', [0.93, 0.69, 0.13], 'CapSize', 8, 'MarkerSize', 6);
errorbar(bin_v_avg4, v_ch4_avg4, v_ch4_se4,...
    'Marker', 'o', 'MarkerFaceColor', [0.85, 0.33, 0.10], 'LineWidth', 1.5, 'Color', [0.85, 0.33, 0.10], 'CapSize', 8, 'MarkerSize', 6);
xlim([0 2]);
ylim([-0.1 0.1]);
xlabel('VPD (kPa)');
ylabel('CO_2-eq flux (mg CO_2-eq m^{-2} s^{-1})');
title('(d) Soil freezing');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);


bn = linspace(0.15, 0.4, 16);
for i = 1:15
    bin_sm1{i} = w.SWC_10_m3m3(w.SWC_10_m3m3 >= bn(i) & w.SWC_10_m3m3 < bn(i+1));
    bin_sm_avg1(i,1) = nanmean(bin_sm1{i});
    sm_ghg1{i} = w.ghg2(w.SWC_10_m3m3 >= bn(i) & w.SWC_10_m3m3 < bn(i+1));
    sm_co2_1{i} = w.co2(w.SWC_10_m3m3 >= bn(i) & w.SWC_10_m3m3 < bn(i+1));
    sm_ch4_1{i} = w.ch4_eq2(w.SWC_10_m3m3 >= bn(i) & w.SWC_10_m3m3 < bn(i+1));
    sm_ghg_avg1(i,1) = nanmean(sm_ghg1{i});
    sm_co2_avg1(i,1) = nanmean(sm_co2_1{i}); 
    sm_ch4_avg1(i,1) = nanmean(sm_ch4_1{i});
    sm_ghg_std1(i,1) = nanstd(sm_ghg1{i});
    sm_co2_std1(i,1) = nanstd(sm_co2_1{i});
    sm_ch4_std1(i,1) = nanstd(sm_ch4_1{i});
    sm_ghg_se1(i,1) = sm_ghg_std1(i,1)/sqrt(length(sm_ghg1{i}));
    sm_co2_se1(i,1) = sm_co2_std1(i,1)/sqrt(length(sm_co2_1{i}));
    sm_ch4_se1(i,1) = sm_ch4_std1(i,1)/sqrt(length(sm_ch4_1{i}));
end

bn = linspace(0.2, 0.52, 16);
for i = 1:15
    bin_sm2{i} = st.SWC_10_m3m3(st.SWC_10_m3m3 >= bn(i) & st.SWC_10_m3m3 < bn(i+1));
    bin_sm_avg2(i,1) = nanmean(bin_sm2{i});
    sm_ghg2{i} = st.ghg2(st.SWC_10_m3m3 >= bn(i) & st.SWC_10_m3m3 < bn(i+1));
    sm_co2_2{i} = st.co2(st.SWC_10_m3m3 >= bn(i) & st.SWC_10_m3m3 < bn(i+1));
    sm_ch4_2{i} = st.ch4_eq2(st.SWC_10_m3m3 >= bn(i) & st.SWC_10_m3m3 < bn(i+1));
    sm_ghg_avg2(i,1) = nanmean(sm_ghg2{i});
    sm_co2_avg2(i,1) = nanmean(sm_co2_2{i});
    sm_ch4_avg2(i,1) = nanmean(sm_ch4_2{i});
    sm_ghg_std2(i,1) = nanstd(sm_ghg2{i});
    sm_co2_std2(i,1) = nanstd(sm_co2_2{i});
    sm_ch4_std2(i,1) = nanstd(sm_ch4_2{i});
    sm_ghg_se2(i,1) = sm_ghg_std2(i,1)/sqrt(length(sm_ghg2{i}));
    sm_co2_se2(i,1) = sm_co2_std2(i,1)/sqrt(length(sm_co2_2{i}));
    sm_ch4_se2(i,1) = sm_ch4_std2(i,1)/sqrt(length(sm_ch4_2{i}));
end

bn = linspace(0.36, 0.54, 16);
for i = 1:15
    bin_sm3{i} = gs.SWC_10_m3m3(gs.SWC_10_m3m3 >= bn(i) & gs.SWC_10_m3m3 < bn(i+1));
    bin_sm_avg3(i,1) = nanmean(bin_sm3{i});
    sm_ghg3{i} = gs.ghg2(gs.SWC_10_m3m3 >= bn(i) & gs.SWC_10_m3m3 < bn(i+1));
    sm_co2_3{i} = gs.co2(gs.SWC_10_m3m3 >= bn(i) & gs.SWC_10_m3m3 < bn(i+1));
    sm_ch4_3{i} = gs.ch4_eq2(gs.SWC_10_m3m3 >= bn(i) & gs.SWC_10_m3m3 < bn(i+1));
    sm_ghg_avg3(i,1) = nanmean(sm_ghg3{i});
    sm_co2_avg3(i,1) = nanmean(sm_co2_3{i});
    sm_ch4_avg3(i,1) = nanmean(sm_ch4_3{i});
    sm_ghg_std3(i,1) = nanstd(sm_ghg3{i});
    sm_co2_std3(i,1) = nanstd(sm_co2_3{i});
    sm_ch4_std3(i,1) = nanstd(sm_ch4_3{i});
    sm_ghg_se3(i,1) = sm_ghg_std3(i,1)/sqrt(length(sm_ghg3{i}));
    sm_co2_se3(i,1) = sm_co2_std3(i,1)/sqrt(length(sm_co2_3{i}));
    sm_ch4_se3(i,1) = sm_ch4_std3(i,1)/sqrt(length(sm_ch4_3{i}));
end

bn = linspace(0.29, 0.54, 16);
for i = 1:15
    bin_sm4{i} = sf.SWC_10_m3m3(sf.SWC_10_m3m3 >= bn(i) & sf.SWC_10_m3m3 < bn(i+1));
    bin_sm_avg4(i,1) = nanmean(bin_sm4{i});
    sm_ghg4{i} = sf.ghg2(sf.SWC_10_m3m3 >= bn(i) & sf.SWC_10_m3m3 < bn(i+1));
    sm_co2_4{i} = sf.co2(sf.SWC_10_m3m3 >= bn(i) & sf.SWC_10_m3m3 < bn(i+1));
    sm_ch4_4{i} = sf.ch4_eq2(sf.SWC_10_m3m3 >= bn(i) & sf.SWC_10_m3m3 < bn(i+1));
    sm_ghg_avg4(i,1) = nanmean(sm_ghg4{i});
    sm_co2_avg4(i,1) = nanmean(sm_co2_4{i});
    sm_ch4_avg4(i,1) = nanmean(sm_ch4_4{i});
    sm_ghg_std4(i,1) = nanstd(sm_ghg4{i});
    sm_co2_std4(i,1) = nanstd(sm_co2_4{i});
    sm_ch4_std4(i,1) = nanstd(sm_ch4_4{i});
    sm_ghg_se4(i,1) = sm_ghg_std4(i,1)/sqrt(length(sm_ghg4{i}));
    sm_co2_se4(i,1) = sm_co2_std4(i,1)/sqrt(length(sm_co2_4{i}));
    sm_ch4_se4(i,1) = sm_ch4_std4(i,1)/sqrt(length(sm_ch4_4{i}));
end


figure(15)
set(gcf, 'Position', [100 100 850 800]);
subplot(2, 2, 1);
hold on; box on;
errorbar(bin_sm_avg1, sm_ghg_avg1, sm_ghg_se1,...
    'Marker', 'o', 'MarkerFaceColor', 'k', 'LineWidth', 1.5, 'Color', 'k', 'CapSize', 8, 'MarkerSize', 6);
errorbar(bin_sm_avg1, sm_co2_avg1, sm_co2_se1,...
    'Marker', 'o', 'MarkerFaceColor', [0.93, 0.69, 0.13], 'LineWidth', 1.5, 'Color', [0.93, 0.69, 0.13], 'CapSize', 8, 'MarkerSize', 6);
errorbar(bin_sm_avg1, sm_ch4_avg1, sm_ch4_se1,...
    'Marker', 'o', 'MarkerFaceColor', [0.85, 0.33, 0.10], 'LineWidth', 1.5, 'Color', [0.85, 0.33, 0.10], 'CapSize', 8, 'MarkerSize', 6);
xlim([0.15 0.4]);
ylim([-0.15 0.1]);
xlabel('SWC 10 cm (m^{3} m^{-3})');
ylabel('CO_2-eq flux (mg CO_2-eq m^{-2} s^{-1})');
legend('net CO_2-eq', 'CO_2', 'CH_4', 'Location', 'southwest', 'FontSize', 11);
title('(a) Winter');
set(gca, 'XTick', -0.05:0.05:0.4, 'FontName', 'Times New Roman', 'FontSize', 14);

subplot(2, 2, 2);
hold on; box on;
errorbar(bin_sm_avg2, sm_ghg_avg2, sm_ghg_se2,...
    'Marker', 'o', 'MarkerFaceColor', 'k', 'LineWidth', 1.5, 'Color', 'k', 'CapSize', 8, 'MarkerSize', 6);
errorbar(bin_sm_avg2, sm_co2_avg2, sm_co2_se2,...
    'Marker', 'o', 'MarkerFaceColor', [0.93, 0.69, 0.13], 'LineWidth', 1.5, 'Color', [0.93, 0.69, 0.13], 'CapSize', 8, 'MarkerSize', 6);
errorbar(bin_sm_avg2, sm_ch4_avg2, sm_ch4_se2, 'Marker', 'o', 'MarkerFaceColor', [0.85, 0.33, 0.10], 'LineWidth', 1.5, 'Color', [0.85, 0.33, 0.10], 'CapSize', 8, 'MarkerSize', 6);
xlim([0.2 0.55]);
ylim([-0.02 0.06]);
xlabel('SWC 10 cm (m^{3} m^{-3})');
ylabel('CO_2-eq flux (mg CO_2-eq m^{-2} s^{-1})');
title('(b) Soil thawing');
set(gca, 'XTick', 0.2:0.05:0.55, 'FontName', 'Times New Roman', 'FontSize', 14);

subplot(2, 2, 3);
hold on; box on;
errorbar(bin_sm_avg3, sm_ghg_avg3, sm_ghg_se3,...
    'Marker', 'o', 'MarkerFaceColor', 'k', 'LineWidth', 1.5, 'Color', 'k', 'CapSize', 8, 'MarkerSize', 6);
errorbar(bin_sm_avg3, sm_co2_avg3, sm_co2_se3,...
    'Marker', 'o', 'MarkerFaceColor', [0.93, 0.69, 0.13], 'LineWidth', 1.5, 'Color', [0.93, 0.69, 0.13], 'CapSize', 8, 'MarkerSize', 6);
errorbar(bin_sm_avg3, sm_ch4_avg3, sm_ch4_se3,...
    'Marker', 'o', 'MarkerFaceColor', [0.85, 0.33, 0.10], 'LineWidth', 1.5, 'Color', [0.85, 0.33, 0.10], 'CapSize', 8, 'MarkerSize', 6);
xlim([0.35 0.55]);
ylim([-0.2 0.1]);
xlabel('SWC 10 cm (m^{3} m^{-3})');
ylabel('CO_2-eq flux (mg CO_2-eq m^{-2} s^{-1})');
title('(c) Growing season');
set(gca,'FontName', 'Times New Roman', 'FontSize', 14);

subplot(2, 2, 4);
hold on; box on;
errorbar(bin_sm_avg4, sm_ghg_avg4, sm_ghg_se4,...
    'Marker', 'o', 'MarkerFaceColor', 'k', 'LineWidth', 1.5, 'Color', 'k', 'CapSize', 8, 'MarkerSize', 6);
errorbar(bin_sm_avg4, sm_co2_avg4, sm_co2_se4,...
    'Marker', 'o', 'MarkerFaceColor', [0.93, 0.69, 0.13], 'LineWidth', 1.5, 'Color', [0.93, 0.69, 0.13], 'CapSize', 8, 'MarkerSize', 6);
errorbar(bin_sm_avg4, sm_ch4_avg4, sm_ch4_se4,...
    'Marker', 'o', 'MarkerFaceColor', [0.85, 0.33, 0.10], 'LineWidth', 1.5, 'Color', [0.85, 0.33, 0.10], 'CapSize', 8, 'MarkerSize', 6);
xlim([0.25 0.55]);
ylim([-0.05 0.1]);
xlabel('SWC 10 cm (m^{3} m^{-3})');
ylabel('CO_2-eq flux (mg CO_2-eq m^{-2} s^{-1})');
title('(d) Soil freezing');
set(gca, 'XTick', 0.25:0.05:0.55, 'FontName', 'Times New Roman', 'FontSize', 14);


bn = linspace(0, 3, 16);
for i = 1:15
    bin_p1{i} = w.PPT_mm(w.PPT_mm >= bn(i) & w.PPT_mm < bn(i+1));
    bin_p_avg1(i,1) = nanmean(bin_p1{i});
    p_ghg1{i} = w.ghg2(w.PPT_mm >= bn(i) & w.PPT_mm < bn(i+1));
    p_co2_1{i} = w.co2(w.PPT_mm >= bn(i) & w.PPT_mm < bn(i+1));
    p_ch4_1{i} = w.ch4_eq2(w.PPT_mm >= bn(i) & w.PPT_mm < bn(i+1));
    p_ghg_avg1(i,1) = nanmean(p_ghg1{i});
    p_co2_avg1(i,1) = nanmean(p_co2_1{i}); 
    p_ch4_avg1(i,1) = nanmean(p_ch4_1{i});
    p_ghg_std1(i,1) = nanstd(p_ghg1{i});
    p_co2_std1(i,1) = nanstd(p_co2_1{i});
    p_ch4_std1(i,1) = nanstd(p_ch4_1{i});
    p_ghg_se1(i,1) = p_ghg_std1(i,1)/sqrt(length(p_ghg1{i}));
    p_co2_se1(i,1) = p_co2_std1(i,1)/sqrt(length(p_co2_1{i}));
    p_ch4_se1(i,1) = p_ch4_std1(i,1)/sqrt(length(p_ch4_1{i}));
end

bn = linspace(0, 5, 16);
for i = 1:15
    bin_p2{i} = st.PPT_mm(st.PPT_mm >= bn(i) & st.PPT_mm < bn(i+1));
    bin_p_avg2(i,1) = nanmean(bin_p2{i});
    p_ghg2{i} = st.ghg2(st.PPT_mm >= bn(i) & st.PPT_mm < bn(i+1));
    p_co2_2{i} = st.co2(st.PPT_mm >= bn(i) & st.PPT_mm < bn(i+1));
    p_ch4_2{i} = st.ch4_eq2(st.PPT_mm >= bn(i) & st.PPT_mm < bn(i+1));
    p_ghg_avg2(i,1) = nanmean(p_ghg2{i});
    p_co2_avg2(i,1) = nanmean(p_co2_2{i});
    p_ch4_avg2(i,1) = nanmean(p_ch4_2{i});
    p_ghg_std2(i,1) = nanstd(p_ghg2{i});
    p_co2_std2(i,1) = nanstd(p_co2_2{i});
    p_ch4_std2(i,1) = nanstd(p_ch4_2{i});
    p_ghg_se2(i,1) = p_ghg_std2(i,1)/sqrt(length(p_ghg2{i}));
    p_co2_se2(i,1) = p_co2_std2(i,1)/sqrt(length(p_co2_2{i}));
    p_ch4_se2(i,1) = p_ch4_std2(i,1)/sqrt(length(p_ch4_2{i}));
end

bn = linspace(0, 18, 16);
for i = 1:15
    bin_p3{i} = gs.PPT_mm(gs.PPT_mm >= bn(i) & gs.PPT_mm < bn(i+1));
    bin_p_avg3(i,1) = nanmean(bin_p3{i});
    p_ghg3{i} = gs.ghg2(gs.PPT_mm >= bn(i) & gs.PPT_mm < bn(i+1));
    p_co2_3{i} = gs.co2(gs.PPT_mm >= bn(i) & gs.PPT_mm < bn(i+1));
    p_ch4_3{i} = gs.ch4_eq2(gs.PPT_mm >= bn(i) & gs.PPT_mm < bn(i+1));
    p_ghg_avg3(i,1) = nanmean(p_ghg3{i});
    p_co2_avg3(i,1) = nanmean(p_co2_3{i});
    p_ch4_avg3(i,1) = nanmean(p_ch4_3{i});
    p_ghg_std3(i,1) = nanstd(p_ghg3{i});
    p_co2_std3(i,1) = nanstd(p_co2_3{i});
    p_ch4_std3(i,1) = nanstd(p_ch4_3{i});
    p_ghg_se3(i,1) = p_ghg_std3(i,1)/sqrt(length(p_ghg3{i}));
    p_co2_se3(i,1) = p_co2_std3(i,1)/sqrt(length(p_co2_3{i}));
    p_ch4_se3(i,1) = p_ch4_std3(i,1)/sqrt(length(p_ch4_3{i}));
end

bn = linspace(0, 3, 16);
for i = 1:15
    bin_p4{i} = sf.PPT_mm(sf.PPT_mm >= bn(i) & sf.PPT_mm < bn(i+1));
    bin_p_avg4(i,1) = nanmean(bin_p4{i});
    p_ghg4{i} = sf.ghg2(sf.PPT_mm >= bn(i) & sf.PPT_mm < bn(i+1));
    p_co2_4{i} = sf.co2(sf.PPT_mm >= bn(i) & sf.PPT_mm < bn(i+1));
    p_ch4_4{i} = sf.ch4_eq2(sf.PPT_mm >= bn(i) & sf.PPT_mm < bn(i+1));
    p_ghg_avg4(i,1) = nanmean(p_ghg4{i});
    p_co2_avg4(i,1) = nanmean(p_co2_4{i});
    p_ch4_avg4(i,1) = nanmean(p_ch4_4{i});
    p_ghg_std4(i,1) = nanstd(p_ghg4{i});
    p_co2_std4(i,1) = nanstd(p_co2_4{i});
    p_ch4_std4(i,1) = nanstd(p_ch4_4{i});
    p_ghg_se4(i,1) = p_ghg_std4(i,1)/sqrt(length(p_ghg4{i}));
    p_co2_se4(i,1) = p_co2_std4(i,1)/sqrt(length(p_co2_4{i}));
    p_ch4_se4(i,1) = p_ch4_std4(i,1)/sqrt(length(p_ch4_4{i}));
end


figure(16)
set(gcf, 'Position', [100 100 850 800]);
subplot(2, 2, 1);
hold on; box on;
errorbar(bin_p_avg1(~isnan(p_ghg_avg1)), p_ghg_avg1(~isnan(p_ghg_avg1)), p_ghg_se1(~isnan(p_ghg_avg1)),...
    'Marker', 'o', 'MarkerFaceColor', 'k', 'LineWidth', 1.5, 'Color', 'k', 'CapSize', 8, 'MarkerSize', 6);
errorbar(bin_p_avg1(~isnan(p_co2_avg1)), p_co2_avg1(~isnan(p_co2_avg1)), p_co2_se1(~isnan(p_co2_avg1)),...
    'Marker', 'o', 'MarkerFaceColor', [0.93, 0.69, 0.13], 'LineWidth', 1.5, 'Color', [0.93, 0.69, 0.13], 'CapSize', 8, 'MarkerSize', 6);
errorbar(bin_p_avg1(~isnan(p_ch4_avg1)), p_ch4_avg1(~isnan(p_ch4_avg1)), p_ch4_se1(~isnan(p_ch4_avg1)),...
    'Marker', 'o', 'MarkerFaceColor', [0.85, 0.33, 0.10], 'LineWidth', 1.5, 'Color', [0.85, 0.33, 0.10], 'CapSize', 8, 'MarkerSize', 6);
xlim([0 3]);
ylim([-0.15 0.1]);
xlabel('PPT (mm)');
ylabel('CO_2-eq flux (mg CO_2-eq m^{-2} s^{-1})');
legend('net CO_2-eq', 'CO_2', 'CH_4', 'FontSize', 11);
title('(a) Winter');
set(gca, 'YTick', -0.15:0.05:0.1, 'FontName', 'Times New Roman', 'FontSize', 14);

subplot(2, 2, 2);
hold on; box on;
errorbar(bin_p_avg2(~isnan(p_ghg_avg2)), p_ghg_avg2(~isnan(p_ghg_avg2)), p_ghg_se2(~isnan(p_ghg_avg2)),...
    'Marker', 'o', 'MarkerFaceColor', 'k', 'LineWidth', 1.5, 'Color', 'k', 'CapSize', 8, 'MarkerSize', 6);
errorbar(bin_p_avg2(~isnan(p_ghg_avg2)), p_co2_avg2(~isnan(p_ghg_avg2)), p_co2_se2(~isnan(p_ghg_avg2)),...
    'Marker', 'o', 'MarkerFaceColor', [0.93, 0.69, 0.13], 'LineWidth', 1.5, 'Color', [0.93, 0.69, 0.13], 'CapSize', 8, 'MarkerSize', 6);
errorbar(bin_p_avg2(~isnan(p_ghg_avg2)), p_ch4_avg2(~isnan(p_ghg_avg2)), p_ch4_se2(~isnan(p_ghg_avg2)),...
    'Marker', 'o', 'MarkerFaceColor', [0.85, 0.33, 0.10], 'LineWidth', 1.5, 'Color', [0.85, 0.33, 0.10], 'CapSize', 8, 'MarkerSize', 6);
xlim([0 5]);
ylim([-0.2 0.2]);
xlabel('PPT (mm)');
ylabel('CO_2-eq flux (mg CO_2-eq m^{-2} s^{-1})');
title('(b) Soil thawing');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);

subplot(2, 2, 3);
hold on; box on;
errorbar(bin_p_avg3(~isnan(p_ghg_avg3)), p_ghg_avg3(~isnan(p_ghg_avg3)), p_ghg_se3(~isnan(p_ghg_avg3)),...
    'Marker', 'o', 'MarkerFaceColor', 'k', 'LineWidth', 1.5, 'Color', 'k', 'CapSize', 8, 'MarkerSize', 6);
errorbar(bin_p_avg3(~isnan(p_ghg_avg3)), p_co2_avg3(~isnan(p_ghg_avg3)), p_co2_se3(~isnan(p_ghg_avg3)),...
    'Marker', 'o', 'MarkerFaceColor', [0.93, 0.69, 0.13], 'LineWidth', 1.5, 'Color', [0.93, 0.69, 0.13], 'CapSize', 8, 'MarkerSize', 6);
errorbar(bin_p_avg3(~isnan(p_ghg_avg3)), p_ch4_avg3(~isnan(p_ghg_avg3)), p_ch4_se3(~isnan(p_ghg_avg3)),...
    'Marker', 'o', 'MarkerFaceColor', [0.85, 0.33, 0.10], 'LineWidth', 1.5, 'Color', [0.85, 0.33, 0.10], 'CapSize', 8, 'MarkerSize', 6);
xlim([0 15]);
ylim([-0.4 0.6]);
xlabel('PPT (mm)');
ylabel('CO_2-eq flux (mg CO_2-eq m^{-2} s^{-1})');
title('(c) Growing season');
set(gca, 'YTick', -0.8:0.2:0.6, 'FontName', 'Times New Roman', 'FontSize', 14);

subplot(2, 2, 4);
hold on; box on;
errorbar(bin_p_avg4(~isnan(p_ghg_avg4)), p_ghg_avg4(~isnan(p_ghg_avg4)), p_ghg_se4(~isnan(p_ghg_avg4)),...
    'Marker', 'o', 'MarkerFaceColor', 'k', 'LineWidth', 1.5, 'Color', 'k', 'CapSize', 8, 'MarkerSize', 6);
errorbar(bin_p_avg4(~isnan(p_ghg_avg4)), p_co2_avg4(~isnan(p_ghg_avg4)), p_co2_se4(~isnan(p_ghg_avg4)),...
    'Marker', 'o', 'MarkerFaceColor', [0.93, 0.69, 0.13], 'LineWidth', 1.5, 'Color', [0.93, 0.69, 0.13], 'CapSize', 8, 'MarkerSize', 6);
errorbar(bin_p_avg4(~isnan(p_ghg_avg4)), p_ch4_avg4(~isnan(p_ghg_avg4)), p_ch4_se4(~isnan(p_ghg_avg4)),...
    'Marker', 'o', 'MarkerFaceColor', [0.85, 0.33, 0.10], 'LineWidth', 1.5, 'Color', [0.85, 0.33, 0.10], 'CapSize', 8, 'MarkerSize', 6);
xlim([0 3]);
ylim([-0.6 0.4]);
xlabel('PPT (mm)');
ylabel('CO_2-eq flux (mg CO_2-eq m^{-2} s^{-1})');
title('(d) Soil freezing');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);


figure(17)
set(gcf, 'Position', [100 100 1200 650]);
subplot(2, 3, 1);
hold on; box on;
errorbar(bin_Rg_avg1, Rg_ghg_avg1, Rg_ghg_se1, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_Rg_avg2, Rg_ghg_avg2, Rg_ghg_se2, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_Rg_avg3, Rg_ghg_avg3, Rg_ghg_se3, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_Rg_avg4, Rg_ghg_avg4, Rg_ghg_se4, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
ylim([-0.4 0.4]);
xlim([0 1500]);
xlabel('R_g (W m^{-2})');
ylabel('net CO_2-eq flux (mg CO_2-eq m^{-2} s^{-1})');
legend('W', 'ST', 'GS', 'SF');
title('(a)');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);

subplot(2, 3, 2);
hold on; box on;
errorbar(bin_Ta_avg1, Ta_ghg_avg1, Ta_ghg_se1, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_Ta_avg2, Ta_ghg_avg2, Ta_ghg_se2, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_Ta_avg3, Ta_ghg_avg3, Ta_ghg_se3, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_Ta_avg4, Ta_ghg_avg4, Ta_ghg_se4, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
ylim([-0.3 0.4]);
xlim([-30 30]);
xlabel('T_{air} (\circC)');
ylabel('net CO_2-eq flux (mg CO_2-eq m^{-2} s^{-1})');
title('(b)');
set(gca, 'YTick', -0.3:0.1:0.4, 'XTick', -30:10:30, 'FontName', 'Times New Roman', 'FontSize', 12);

subplot(2, 3, 3);
hold on; box on;
errorbar(bin_p_avg1(~isnan(p_ghg_avg1)), p_ghg_avg1(~isnan(p_ghg_avg1)), p_ghg_se1(~isnan(p_ghg_avg1)),...
    'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_p_avg2(~isnan(p_ghg_avg2)), p_ghg_avg2(~isnan(p_ghg_avg2)), p_ghg_se2(~isnan(p_ghg_avg2)),...
    'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_p_avg3(~isnan(p_ghg_avg3)), p_ghg_avg3(~isnan(p_ghg_avg3)), p_ghg_se3(~isnan(p_ghg_avg3)),...
    'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_p_avg4(~isnan(p_ghg_avg4)), p_ghg_avg4(~isnan(p_ghg_avg4)), p_ghg_se4(~isnan(p_ghg_avg4)),...
    'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
ylim([-0.2 0.6]);
xlabel('PPT (mm)');
ylabel('net CO_2-eq flux (mg CO_2-eq m^{-2} s^{-1})');
title('(c)');
set(gca, 'YTick', -0.6:0.2:0.6, 'FontName', 'Times New Roman', 'FontSize', 12);
% breakxaxis([9 13]);

subplot(2, 3, 4);
hold on; box on;
errorbar(bin_v_avg1, v_ghg_avg1, v_ghg_se1, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_v_avg2, v_ghg_avg2, v_ghg_se2, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_v_avg3, v_ghg_avg3, v_ghg_se3, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_v_avg4, v_ghg_avg4, v_ghg_se4, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
ylim([-0.3 0.2]);
xlim([0 2.5]);
xlabel('VPD (kPa)');
ylabel('net CO_2-eq flux (mg CO_2-eq m^{-2} s^{-1})');
title('(d)');
set(gca, 'XTick', 0:0.5:2.5, 'FontName', 'Times New Roman', 'FontSize', 12);

subplot(2, 3, 5);
hold on; box on;
errorbar(bin_Ts_avg1, Ts_ghg_avg1, Ts_ghg_se1, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_Ts_avg2, Ts_ghg_avg2, Ts_ghg_se2, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_Ts_avg3, Ts_ghg_avg3, Ts_ghg_se3, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_Ts_avg4, Ts_ghg_avg4, Ts_ghg_se4, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
ylim([-0.2 0.2]);
xlim([-5 25]);
xlabel('T_{soil} 10 cm (\circC)');
ylabel('net CO_2-eq flux (mg CO_2-eq m^{-2} s^{-1})');
title('(e)');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);

subplot(2, 3, 6);
hold on; box on;
errorbar(bin_sm_avg1, sm_ghg_avg1, sm_ghg_se1, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_sm_avg2, sm_ghg_avg2, sm_ghg_se2, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_sm_avg3, sm_ghg_avg3, sm_ghg_se3, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_sm_avg4, sm_ghg_avg4, sm_ghg_se4, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
ylim([-0.2 0.2]);
xlim([0.1 0.6]);
xlabel('SWC 10 cm (m^{3} m^{-3})');
ylabel('net CO_2-eq flux (mg CO_2-eq m^{-2} s^{-1})');
title('(f)');
set(gca, 'XTick', 0.1:0.1:0.6, 'FontName', 'Times New Roman', 'FontSize', 12);


figure(18)
set(gcf, 'Position', [100 100 1200 650]);
subplot(2, 3, 1);
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

subplot(2, 3, 2);
hold on; box on;
errorbar(bin_Ta_avg1, Ta_co2_avg1, Ta_co2_se1, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_Ta_avg2, Ta_co2_avg2, Ta_co2_se2, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_Ta_avg3, Ta_co2_avg3, Ta_co2_se3, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_Ta_avg4, Ta_co2_avg4, Ta_co2_se4, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
ylim([-0.5 0.5]);
xlim([-30 30]);
xlabel('T_{air} (\circC)');
ylabel('CO_2 flux (mg CO_2 m^{-2} s^{-1})');
title('(b)');
set(gca, 'XTick', -30:10:30, 'FontName', 'Times New Roman', 'FontSize', 12);

subplot(2, 3, 3);
hold on; box on;
errorbar(bin_p_avg1(~isnan(p_co2_avg1)), p_co2_avg1(~isnan(p_co2_avg1)), p_co2_se1(~isnan(p_co2_avg1)),...
    'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_p_avg2(~isnan(p_co2_avg2)), p_co2_avg2(~isnan(p_co2_avg2)), p_ghg_se2(~isnan(p_co2_avg2)),...
    'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_p_avg3(~isnan(p_co2_avg3)), p_co2_avg3(~isnan(p_co2_avg3)), p_ghg_se3(~isnan(p_co2_avg3)),...
    'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_p_avg4(~isnan(p_co2_avg4)), p_co2_avg4(~isnan(p_co2_avg4)), p_ghg_se4(~isnan(p_co2_avg4)),...
    'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
ylim([-1 0.5]);
xlabel('PPT (mm)');
ylabel('CO_2 flux (mg CO_2 m^{-2} s^{-1})');
title('(c)');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);
breakxaxis([9 13]);

subplot(2, 3, 4);
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

subplot(2, 3, 5);
hold on; box on;
errorbar(bin_Ts_avg1, Ts_co2_avg1, Ts_co2_se1, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_Ts_avg2, Ts_co2_avg2, Ts_co2_se2, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_Ts_avg3, Ts_co2_avg3, Ts_co2_se3, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_Ts_avg4, Ts_co2_avg4, Ts_co2_se4, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
ylim([-0.3 0.1]);
xlim([-5 25]);
xlabel('T_{soil} 10 cm (\circC)');
ylabel('CO_2 flux (mg CO_2 m^{-2} s^{-1})');
title('(e)');
set(gca, 'YTick', -0.3:0.1:0.1, 'FontName', 'Times New Roman', 'FontSize', 12);

subplot(2, 3, 6);
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


figure(19)
set(gcf, 'Position', [100 100 1200 650]);
subplot(2, 3, 1);
hold on; box on;
errorbar(bin_Rg_avg1, Rg_ch4_avg1, Rg_ch4_se1, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_Rg_avg2, Rg_ch4_avg2, Rg_ch4_se2, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_Rg_avg3, Rg_ch4_avg3, Rg_ch4_se3, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_Rg_avg4, Rg_ch4_avg4, Rg_ch4_se4, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
ylim([0 0.1]);
xlim([0 1500]);
xlabel('R_g (W m^{-2})');
ylabel('CH_4 flux (mg CO_2-eq m^{-2} s^{-1})');
legend('W', 'ST', 'GS', 'SF', 'Location', 'northeast', 'Orientation', 'horizontal', 'FontSize' ,10, 'NumColumns', 2);
title('(a)');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);

subplot(2, 3, 2);
hold on; box on;
errorbar(bin_Ta_avg1, Ta_ch4_avg1, Ta_ch4_se1, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_Ta_avg2, Ta_ch4_avg2, Ta_ch4_se2, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_Ta_avg3, Ta_ch4_avg3, Ta_ch4_se3, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_Ta_avg4, Ta_ch4_avg4, Ta_ch4_se4, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
ylim([-0.05 0.15]);
xlim([-30 30]);
xlabel('T_{air} (\circC)');
ylabel('CH_4 flux (mg CO_2-eq m^{-2} s^{-1})');
title('(b)');
set(gca, 'XTick', -30:10:30, 'FontName', 'Times New Roman', 'FontSize', 12);

subplot(2, 3, 3);
hold on; box on;
errorbar(bin_p_avg1(~isnan(p_ch4_avg1)), p_ch4_avg1(~isnan(p_ch4_avg1)), p_ch4_se1(~isnan(p_ch4_avg1)),...
    'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_p_avg2(~isnan(p_ch4_avg2)), p_ch4_avg2(~isnan(p_ch4_avg2)), p_ghg_se2(~isnan(p_ch4_avg2)),...
    'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_p_avg3(~isnan(p_ch4_avg3)), p_ch4_avg3(~isnan(p_ch4_avg3)), p_ghg_se3(~isnan(p_ch4_avg3)),...
    'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_p_avg4(~isnan(p_ch4_avg4)), p_ch4_avg4(~isnan(p_ch4_avg4)), p_ghg_se4(~isnan(p_ch4_avg4)),...
    'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
ylim([-0.05 0.2]);
xlabel('PPT (mm)');
ylabel('CH_4 flux (mg CO_2-eq m^{-2} s^{-1})');
title('(c)');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);
% breakxaxis([9 13]);

subplot(2, 3, 4);
hold on; box on;
errorbar(bin_v_avg1, v_ch4_avg1, v_ch4_se1, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_v_avg2, v_ch4_avg2, v_ch4_se2, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_v_avg3, v_ch4_avg3, v_ch4_se3, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_v_avg4, v_ch4_avg4, v_ch4_se4, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
ylim([0 0.1]);
xlim([0 2.5]);
xlabel('VPD (kPa)');
ylabel('CH_4 flux (mg CO_2-eq m^{-2} s^{-1})');
title('(d)');
set(gca, 'XTick', 0:0.5:2.5, 'FontName', 'Times New Roman', 'FontSize', 12);

subplot(2, 3, 5);
hold on; box on;
errorbar(bin_Ts_avg1, Ts_ch4_avg1, Ts_ch4_se1, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_Ts_avg2, Ts_ch4_avg2, Ts_ch4_se2, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_Ts_avg3, Ts_ch4_avg3, Ts_ch4_se3, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_Ts_avg4, Ts_ch4_avg4, Ts_ch4_se4, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
ylim([-0.05 0.15]);
xlim([-5 25]);
xlabel('T_{soil} 10 cm (\circC)');
ylabel('CH_4 flux (mg CO_2-eq m^{-2} s^{-1})');
title('(e)');
set(gca, 'YTick', -0.05:0.05:0.15, 'FontName', 'Times New Roman', 'FontSize', 12);

subplot(2, 3, 6);
hold on; box on;
errorbar(bin_sm_avg1, sm_ch4_avg1, sm_ch4_se1, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_sm_avg2, sm_ch4_avg2, sm_ch4_se2, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_sm_avg3, sm_ch4_avg3, sm_ch4_se3, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
errorbar(bin_sm_avg4, sm_ch4_avg4, sm_ch4_se4, 'Marker', 'o', 'CapSize', 6, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');
ylim([0 0.1]);
xlim([0.1 0.6]);
xlabel('SWC 10 cm (m^{3} m^{-3})');
ylabel('CH_4 flux (mg CO_2-eq m^{-2} s^{-1})');
title('(f)');
set(gca, 'XTick', 0.1:0.1:0.6, 'FontName', 'Times New Roman', 'FontSize', 12);



x3 = find(ismember(mydat.TIMESTAMP, '10/02/2014 00:30')); %update the beginning time
y3 = find(ismember(mydat.TIMESTAMP, '12/03/2014 00:00')); %update the ending time
mydat3 = mydat(x3:y3, :);

rg_3 = reshape(w1.Rg_Wm2, 48, []);
rg_3_avg = nanmean(rg_3');
rg_3_std = nanstd(rg_3');
rg_3_se = rg_3_std/length(rg_3);

ta_3 = reshape(w1.Ta_degC, 48, []);
ta_3_avg = nanmean(ta_3');
ta_3_std = nanstd(ta_3');
ta_3_se = ta_3_std/length(ta_3);

v_3 = reshape(w1.VPD_kPa, 48, []);
v_3_avg = nanmean(v_3');
v_3_std = nanstd(v_3');
v_3_se = v_3_std/length(v_3);

ts_3 = reshape(w1.Soil_10_degC, 48, []);
ts_3_avg = nanmean(ts_3');
ts_3_std = nanstd(ts_3');
ts_3_se = ts_3_std/length(ts_3);

sm_3 = reshape(w1.SWC_10_m3m3,48,[]);
sm_3_avg = nanmean(sm_3');
sm_3_std = nanstd(sm_3');
sm_3_se = sm_3_std/length(sm_3);

co2_3 = reshape(w1.co2_flux_filled_umolm2s*44/1000, 48, []);
co2_3_avg = nanmean(co2_3');
co2_3_std = nanstd(co2_3');
co2_3_se = co2_3_std/length(co2_3);

ch4_3 = reshape(w1.ch4_flux_filled_umolm2s*16*28/1000, 48, []);
ch4_3_avg = nanmean(ch4_3');
ch4_3_std = nanstd(ch4_3');
ch4_3_se = ch4_3_std/length(ch4_3);

ghg_3 = reshape((w1.co2_flux_filled_umolm2s*44/1000+w1.ch4_flux_filled_umolm2s*16*28/1000), 48, []);
ghg_3_avg = nanmean(ghg_3');
ghg_3_std = nanstd(ghg_3');
ghg_3_se = ghg_3_std/length(ghg_3);


x4 = find(ismember(mydat.TIMESTAMP, '10/08/2015 00:30')); %update the beginning time
y4 = find(ismember(mydat.TIMESTAMP, '12/17/2015 00:00')); %update the ending time
mydat4 = mydat(x4:y4, :);

rg_4 = reshape(w2.Rg_Wm2, 48, []);
rg_4_avg = nanmean(rg_4');
rg_4_std = nanstd(rg_4');
rg_4_se = rg_4_std/length(rg_4);

ta_4 = reshape(w2.Ta_degC, 48, []);
ta_4_avg = nanmean(ta_4');
ta_4_std = nanstd(ta_4');
ta_4_se = ta_4_std/length(ta_4);

v_4 = reshape(w2.VPD_kPa, 48, []);
v_4_avg = nanmean(v_4');
v_4_std = nanstd(v_4');
v_4_se = v_4_std/length(v_4);

ts_4 = reshape(w2.Soil_10_degC, 48, []);
ts_4_avg = nanmean(ts_4');
ts_4_std = nanstd(ts_4');
ts_4_se = ts_4_std/length(ts_4);

sm_4 = reshape(w2.SWC_10_m3m3, 48, []);
sm_4_avg = nanmean(sm_4');
sm_4_std = nanstd(sm_4');
sm_4_se = sm_4_std/length(sm_4);

co2_4 = reshape(w2.co2_flux_filled_umolm2s*44/1000, 48, []);
co2_4_avg = nanmean(co2_4');
co2_4_std = nanstd(co2_4');
co2_4_se = co2_4_std/length(co2_4);

ch4_4 = reshape(w2.ch4_flux_filled_umolm2s*16*28/1000, 48, []);
ch4_4_avg = nanmean(ch4_4');
ch4_4_std = nanstd(ch4_4');
ch4_4_se = ch4_4_std/length(ch4_4);

ghg_4 = reshape((w2.co2_flux_filled_umolm2s*44/1000+w2.ch4_flux_filled_umolm2s*16*28/1000), 48, []);
ghg_4_avg = nanmean(ghg_4');
ghg_4_std = nanstd(ghg_4');
ghg_4_se = ghg_4_std/length(ghg_4);


x5 = find(ismember(mydat.TIMESTAMP, '05/01/2016 00:30')); %update the beginning time
y5 = find(ismember(mydat.TIMESTAMP, '08/01/2016 00:00')); %update the ending time
mydat5 = mydat(x5:y5, :);

rg_5 = reshape(mydat5.Rg_Wm2, 48 ,[]);
rg_5_avg = nanmean(rg_5');
rg_5_std = nanstd(rg_5');
rg_5_se = rg_5_std/length(rg_5);

ta_5 = reshape(mydat5.Ta_degC, 48, []);
ta_5_avg = nanmean(ta_5');
ta_5_std = nanstd(ta_5');
ta_5_se = ta_5_std/length(ta_5);

v_5 = reshape(mydat5.VPD_kPa, 48, []);
v_5_avg = nanmean(v_5');
v_5_std = nanstd(v_5');
v_5_se = v_5_std/length(v_5);

ts_5 = reshape(mydat5.Soil_10_degC, 48, []);
ts_5_avg = nanmean(ts_5');
ts_5_std = nanstd(ts_5');
ts_5_se = ts_5_std/length(ts_5);

sm_5 = reshape(mydat5.SWC_10_m3m3, 48, []);
sm_5_avg = nanmean(sm_5');
sm_5_std = nanstd(sm_5');
sm_5_se = sm_5_std/length(sm_5);

co2_5 = reshape(mydat5.co2_flux_filled_umolm2s*44/1000, 48, []);
co2_5_avg = nanmean(co2_5');
co2_5_std = nanstd(co2_5');
co2_5_se = co2_5_std/length(co2_5);

ch4_5 = reshape(mydat5.ch4_flux_filled_umolm2s*16*28/1000, 48, []);
ch4_5_avg = nanmean(ch4_5');
ch4_5_std = nanstd(ch4_5');
ch4_5_se = ch4_5_std/length(ch4_5);

ghg_5 = reshape((mydat5.co2_flux_filled_umolm2s*44/1000+mydat5.ch4_flux_filled_umolm2s*16*28/1000), 48, []);
ghg_5_avg = nanmean(ghg_5');
ghg_5_std = nanstd(ghg_5');
ghg_5_se = ghg_5_std/length(ghg_5);


figure(20)
set(gcf, 'Position' ,[20 20 1200 1000]);
subplot(3, 3, 1);
hold on; box on;
shadedErrorBar(1:48, ghg_3_avg, ghg_3_se, 'lineProps', '-r');
h1 = plot(1:48, ghg_3_avg, 'r', 'LineWidth', 1.5);
shadedErrorBar(1:48, ghg_4_avg, ghg_4_se, 'lineProps', '-b');
h2 = plot(1:48, ghg_4_avg, 'b', 'LineWidth', 1.5);
shadedErrorBar(1:48, ghg_5_avg, ghg_5_se, 'lineProps', '-g');
h3 = plot(1:48, ghg_5_avg, 'g', 'LineWidth', 1.5);
yline(0, '--k');
legend([h1 h2 h3], {'2014','2015','2016'}, 'Location', 'southeast', 'FontSize',9);
xlim([0 48]);
% ylim([-0.1 0.2])
xlabel('Hour');
ylabel('net CO_2-eq flux (mg CO_2-eq m^{-2} s^{-1})');
title('(a)');
set(gca, 'XTick', 0:12:48, 'XTickLabel', {'0','6','12','18','24'}, 'FontName', 'Times New Roman', 'FontSize', 14);

subplot(3, 3, 2);
hold on; box on;
shadedErrorBar(1:48, co2_3_avg, co2_3_se, 'lineProps', '-r');
h1 = plot(1:48, co2_3_avg, 'r', 'LineWidth', 1.5);
shadedErrorBar(1:48, co2_4_avg, co2_4_se, 'lineProps', '-b');
h2 = plot(1:48, co2_4_avg, 'b', 'LineWidth', 1.5);
shadedErrorBar(1:48, co2_5_avg, co2_5_se, 'lineProps', '-g');
h3 = plot(1:48, co2_5_avg, 'g', 'LineWidth', 1.5);
yline(0, '--k');
legend([h1 h2 h3], {'2014','2015','2016'}, 'Location', 'southeast', 'FontSize', 9);
xlim([0 48]);
% ylim([-0.1 0.2])
xlabel('Hour');
ylabel('CO_2 flux (mg CO_2 m^{-2} s^{-1})');
title('(b)');
set(gca, 'XTick', 0:12:48, 'XTickLabel', {'0','6','12','18','24'}, 'FontName', 'Times New Roman', 'FontSize', 14);

subplot(3, 3, 3);
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
% ylim([0 0.02])
xlabel('Hour');
ylabel('CH_4 flux (mg CO_2-eq m^{-2} s^{-1})');
title('(c)');
set(gca,'XTick', 0:12:48, 'XTickLabel', {'0','6','12','18','24'}, 'FontName', 'Times New Roman', 'FontSize', 14);

subplot(3, 3, 4);
hold on; box on;
shadedErrorBar(1:48, rg_3_avg, rg_3_se, 'lineProps', '-r');
h1 = plot(1:48, rg_3_avg, 'r', 'LineWidth', 1.5);
shadedErrorBar(1:48, rg_4_avg, rg_4_se, 'lineProps', '-b');
h2 = plot(1:48, rg_4_avg, 'b', 'LineWidth', 1.5);
shadedErrorBar(1:48, rg_5_avg, rg_5_se, 'lineProps', '-g');
h3 = plot(1:48, rg_5_avg, 'g', 'LineWidth', 1.5);
legend([h1 h2 h3], {'2014','2015','2016'}, 'Location', 'northwest', 'FontSize', 9);
xlim([0 48]);
% ylim([0 800]);
xlabel('Hour');
ylabel('R_g (W m^{-2})');
title('(d)');
set(gca, 'XTick', 0:12:48, 'XTickLabel', {'0','6','12','18','24'}, 'FontName', 'Times New Roman', 'FontSize', 14);

subplot(3, 3, 5);
hold on; box on;
shadedErrorBar(1:48, ta_3_avg, ta_3_se, 'lineProps', '-r');
h1 = plot(1:48, ta_3_avg, 'r', 'LineWidth', 1.5);
shadedErrorBar(1:48, ta_4_avg, ta_4_se, 'lineProps', '-b');
h2 = plot(1:48, ta_4_avg, 'b', 'LineWidth', 1.5);
shadedErrorBar(1:48, ta_5_avg, ta_5_se, 'lineProps', '-g');
h3 = plot(1:48, ta_5_avg, 'g', 'LineWidth', 1.5);
legend([h1 h2 h3], {'2014','2015','2016'}, 'Location', 'northwest', 'FontSize', 9);
xlim([0 48]);
% ylim([2 16]);
xlabel('Hour');
ylabel('T_{air} (\circC)');
title('(e)');
set(gca, 'XTick', 0:12:48, 'XTickLabel', {'0','6','12','18','24'}, 'FontName', 'Times New Roman', 'FontSize', 14);

subplot(3, 3, 6);
hold on; box on;
h = bar(2014:2016, [nansum(mydat3.PPT_mm), nansum(mydat4.PPT_mm), nansum(mydat5.PPT_mm)], 'BarWidth', 0.4);
h.FaceColor = 'flat';
h.CData(1,:) = [1 0 0];
h.CData(2,:) = [0 0 1]';
h.CData(3,:) = [0 1 0];
ylabel('PPT (mm)');
title('(f)');
set(gca, 'XTick', 2014:2016, 'FontName', 'Times New Roman', 'FontSize', 14);

subplot(3, 3, 7);
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

subplot(3, 3, 8);
hold on; box on;
shadedErrorBar(1:48, ts_3_avg, ts_3_se, 'lineProps', '-r');
h1 = plot(1:48, ts_3_avg, 'r', 'LineWidth', 1.5);
shadedErrorBar(1:48, ts_4_avg, ts_4_se, 'lineProps', '-b');
h2 = plot(1:48, ts_4_avg, 'b', 'LineWidth', 1.5);
shadedErrorBar(1:48, ts_5_avg, ts_5_se, 'lineProps', '-g');
h3 = plot(1:48, ts_5_avg, 'g', 'LineWidth', 1.5);
legend([h1 h2 h3], {'2014','2015','2016'}, 'Location', 'southeast', 'FontSize', 9);
xlim([0 48]);
% ylim([8 16]);
xlabel('Hour');
ylabel('T_{soil} 10cm (\circC)');
title('(h)');
set(gca, 'XTick', 0:12:48, 'XTickLabel', {'0','6','12','18','24'}, 'FontName', 'Times New Roman', 'FontSize', 14);

subplot(3, 3, 9);
hold on; box on;
shadedErrorBar(1:48, sm_3_avg, sm_3_se, 'lineProps', '-r');
h1 = plot(1:48, sm_3_avg, 'r', 'LineWidth', 1.5);
shadedErrorBar(1:48, sm_4_avg, sm_4_se, 'lineProps', '-b');
h2 = plot(1:48, sm_4_avg, 'b', 'LineWidth', 1.5);
shadedErrorBar(1:48, sm_5_avg, sm_5_se, 'lineProps', '-g');
h3 = plot(1:48, sm_5_avg, 'g', 'LineWidth', 1.5);
legend([h1 h2 h3], {'2014','2015','2016'}, 'Location', 'east', 'FontSize', 9);
xlim([0 48]);
% ylim([0.45 0.5]);
xlabel('Hour');
ylabel('SWC 10cm (m^3 m^{-3})');
title('(i)');
set(gca, 'XTick', 0:12:48, 'XTickLabel', {'0','6','12','18','24'}, 'FontName', 'Times New Roman', 'FontSize', 14);


figure(21)
set(gcf, 'Position' ,[20 20 1200 1000]);
subplot(3, 3, 1);
hold on; box on;
shadedErrorBar(1:48, ghg_3_avg, ghg_3_se, 'lineProps', '-r');
h1 = plot(1:48, ghg_3_avg, 'r', 'LineWidth', 1.5);
shadedErrorBar(1:48, ghg_4_avg, ghg_4_se, 'lineProps', '-b');
h2 = plot(1:48, ghg_4_avg, 'b', 'LineWidth', 1.5);
yline(0, '--k');
legend([h1 h2], {'2014','2015'}, 'Location', 'southeast', 'FontSize',9);
xlim([0 48]);
% ylim([-0.1 0.2])
xlabel('Hour');
ylabel('net CO_2-eq flux (mg CO_2-eq m^{-2} s^{-1})');
title('(a)');
set(gca, 'XTick', 0:12:48, 'XTickLabel', {'0','6','12','18','24'}, 'FontName', 'Times New Roman', 'FontSize', 14);

subplot(3, 3, 2);
hold on; box on;
shadedErrorBar(1:48, co2_3_avg, co2_3_se, 'lineProps', '-r');
h1 = plot(1:48, co2_3_avg, 'r', 'LineWidth', 1.5);
shadedErrorBar(1:48, co2_4_avg, co2_4_se, 'lineProps', '-b');
h2 = plot(1:48, co2_4_avg, 'b', 'LineWidth', 1.5);
yline(0, '--k');
legend([h1 h2], {'2014','2015'}, 'Location', 'southeast', 'FontSize', 9);
xlim([0 48]);
% ylim([-0.1 0.2])
xlabel('Hour');
ylabel('CO_2 flux (mg CO_2 m^{-2} s^{-1})');
title('(b)');
set(gca, 'XTick', 0:12:48, 'XTickLabel', {'0','6','12','18','24'}, 'FontName', 'Times New Roman', 'FontSize', 14);

subplot(3, 3, 3);
hold on; box on;
shadedErrorBar(1:48, ch4_3_avg, ch4_3_se, 'lineProps', '-r');
h1 = plot(1:48, ch4_3_avg, 'r', 'LineWidth', 1.5);
shadedErrorBar(1:48, ch4_4_avg, ch4_4_se, 'lineProps', '-b');
h2 = plot(1:48, ch4_4_avg, 'b', 'LineWidth', 1.5);
yline(0, '--k');
legend([h1 h2], {'2014','2015'}, 'Location', 'southeast', 'FontSize', 9);
xlim([0 48]);
% ylim([0 0.02])
xlabel('Hour');
ylabel('CH_4 flux (mg CO_2-eq m^{-2} s^{-1})');
title('(c)');
set(gca,'XTick', 0:12:48, 'XTickLabel', {'0','6','12','18','24'}, 'FontName', 'Times New Roman', 'FontSize', 14);

subplot(3, 3, 4);
hold on; box on;
shadedErrorBar(1:48, rg_3_avg, rg_3_se, 'lineProps', '-r');
h1 = plot(1:48, rg_3_avg, 'r', 'LineWidth', 1.5);
shadedErrorBar(1:48, rg_4_avg, rg_4_se, 'lineProps', '-b');
h2 = plot(1:48, rg_4_avg, 'b', 'LineWidth', 1.5);
legend([h1 h2], {'2014','2015'}, 'Location', 'northwest', 'FontSize', 9);
xlim([0 48]);
% ylim([0 800]);
xlabel('Hour');
ylabel('R_g (W m^{-2})');
title('(d)');
set(gca, 'XTick', 0:12:48, 'XTickLabel', {'0','6','12','18','24'}, 'FontName', 'Times New Roman', 'FontSize', 14);

subplot(3, 3, 5);
hold on; box on;
shadedErrorBar(1:48, ta_3_avg, ta_3_se, 'lineProps', '-r');
h1 = plot(1:48, ta_3_avg, 'r', 'LineWidth', 1.5);
shadedErrorBar(1:48, ta_4_avg, ta_4_se, 'lineProps', '-b');
h2 = plot(1:48, ta_4_avg, 'b', 'LineWidth', 1.5);
legend([h1 h2 h3], {'2014','2015'}, 'Location', 'northwest', 'FontSize', 9);
xlim([0 48]);
% ylim([2 16]);
xlabel('Hour');
ylabel('T_{air} (\circC)');
title('(e)');
set(gca, 'XTick', 0:12:48, 'XTickLabel', {'0','6','12','18','24'}, 'FontName', 'Times New Roman', 'FontSize', 14);

subplot(3, 3, 6);
hold on; box on;
h = bar(2014:2015, [nansum(w1.PPT_mm), nansum(w2.PPT_mm)], 'BarWidth', 0.4);
h.FaceColor = 'flat';
h.CData(1,:) = [1 0 0];
h.CData(2,:) = [0 0 1]';
ylabel('PPT (mm)');
title('(f)');
set(gca, 'XTick', 2014:2015, 'FontName', 'Times New Roman', 'FontSize', 14);

subplot(3, 3, 7);
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

subplot(3, 3, 8);
hold on; box on;
shadedErrorBar(1:48, ts_3_avg, ts_3_se, 'lineProps', '-r');
h1 = plot(1:48, ts_3_avg, 'r', 'LineWidth', 1.5);
shadedErrorBar(1:48, ts_4_avg, ts_4_se, 'lineProps', '-b');
h2 = plot(1:48, ts_4_avg, 'b', 'LineWidth', 1.5);
legend([h1 h2], {'2014','2015'}, 'Location', 'southeast', 'FontSize', 9);
xlim([0 48]);
% ylim([8 16]);
xlabel('Hour');
ylabel('T_{soil} 10cm (\circC)');
title('(h)');
set(gca, 'XTick', 0:12:48, 'XTickLabel', {'0','6','12','18','24'}, 'FontName', 'Times New Roman', 'FontSize', 14);

subplot(3, 3, 9);
hold on; box on;
shadedErrorBar(1:48, sm_3_avg, sm_3_se, 'lineProps', '-r');
h1 = plot(1:48, sm_3_avg, 'r', 'LineWidth', 1.5);
shadedErrorBar(1:48, sm_4_avg, sm_4_se, 'lineProps', '-b');
h2 = plot(1:48, sm_4_avg, 'b', 'LineWidth', 1.5);
legend([h1 h2], {'2014','2015'}, 'Location', 'east', 'FontSize', 9);
xlim([0 48]);
% ylim([0.45 0.5]);
xlabel('Hour');
ylabel('SWC 10cm (m^3 m^{-3})');
title('(i)');
set(gca, 'XTick', 0:12:48, 'XTickLabel', {'0','6','12','18','24'}, 'FontName', 'Times New Roman', 'FontSize', 14);