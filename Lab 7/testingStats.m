% Code for PHY 105N Lab 7 Data Analysis
% This script reads q/m data from an Excel file, creates a scatterplot,

% === 1. Load data ===
T = readtable("PHY105N_Lab7Data.xlsx");
%TableProperties = T.Properties.VariableNames;

% Change this to the actual variable name in the table:
qm = T.q_m;   % e.g. T.q_m, T.q_m_, T.q_m_C_kg_, etc.

n = numel(qm);

% === 2. Scatterplot (q/m vs trial number) ===
trials = (1:n)';

figure;
scatter(trials, qm, 'filled');
xlabel('Trial number');
ylabel('q/m (C/kg)');
title('q/m measurements');
grid on;

% === 3. t-score vs accepted value ===
qm_mean = mean(qm);
qm_std  = std(qm, 0);         % 0 => sample std (divide by n-1)
se_mean = qm_std / sqrt(n);

qm_accepted = 1.76e11;        % accepted q/m for electron (C/kg)

t_score = (qm_mean - qm_accepted) / se_mean;

fprintf('=======Data Analysis Results=======\n');
fprintf('Mean q/m: %g\n', qm_mean);
fprintf('Std dev: %g\n', qm_std);
fprintf('Standard error: %g\n', se_mean);
fprintf('t-score vs accepted: %g\n', t_score);

% === 4. Chi-square vs accepted value ===
% Option A: assume same sigma for each point = sample std:
sigma = qm_std;

chi2 = sum((qm - qm_accepted).^2 / sigma^2);
dof = n - 1;
chi2_red = chi2 / dof;

fprintf('Chi-square: %g\n', chi2);
fprintf('Reduced chi-square: %g\n', chi2_red);
