% T-score testing against the data
% tscoreTesting.m
% This script performs t-score testing on capacitor data from an Excel file
% to compare two sets of measurements: Capacitance vs Plate Area and
% Capacitance vs Plate Separation.
%
% PHY 105N, Abdon Morales, am226923
T = readtable('CapacitorFinalProjectDataSet.xlsx', 'VariableNamingRule', 'preserve');

%Data
C_sep = T.("Capacitance (pF)")(7:12); % Capacitance vs Plate Separation
C_area = T.("Capacitance (pF)")(1:6); % Capacitance vs Plate Area

% Hypothesized mean
mu0 = 0;

% Perform t-test for each set
n_sep = numel(C_sep);
mean_sep = mean(C_sep);
sd_sep = std(C_sep, 1); % The standard deviation for population, for sample use std(C_sep, 0)
t_sep = (mean_sep - mu0) / (sd_sep / sqrt(n_sep));
df_sep = n_sep - 1;
p_value_sep = 2 * (1 - tcdf(abs(t_sep), df_sep));

% Dump results
fprintf('Seperation: mean=%.4f, sd=%.4f, n=%d, t=%.4f, p=%.4f\n', ...
    mean_sep, sd_sep, n_sep, t_sep, p_value_sep);

% Area dataset
n_area = numel(C_area);
mean_area = mean(C_area);
sd_area = std(C_area, 1); % The standard deviation for population, for sample use std(C_area, 0)
t_area = (mean_area - mu0) / (sd_area / sqrt(n_area));
df_area = n_area - 1;
p_value_area = 2 * (1 - tcdf(abs(t_area), df_area));
% Dump results
fprintf('Area: mean=%.4f, sd=%.4f, n=%d, t=%.4f, p=%.4f\n', ...
    mean_area, sd_area, n_area, t_area, p_value_area);

% Two sample t-test (Independent)
[~, p_value_two_sample, ci, stats] = ttest2(C_sep, C_area);
fprintf('Two-sample t-test: t=%.4f, df=%d, p=%.4f\n', ...
    stats.tstat, stats.df, p_value_two_sample);
