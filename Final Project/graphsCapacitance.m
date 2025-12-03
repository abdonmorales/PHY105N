% graphsCapacitance.m
% This script reads capacitor data from an Excel file and generates
% two plots: Capacitance vs Plate Area and Capacitance vs Plate Separation.
%
% PHY 105N, Abdon Morales, am226923
% Load data from Excel
data = readtable('CapacitorFinalProjectDataSet.xlsx', 'VariableNamingRule', 'preserve');

% Extract subsets
area_trials = data(1:6, :);       % Rows 2–7
sep_trials = data(7:12, :);       % Rows 8–13

% Plot Capacitance vs Plate Area
figure;
plot(area_trials.("Plate Area (mm^2)"), area_trials.("Capacitance (pF)"), 'o-', 'LineWidth', 2);
xlabel('Plate Area (mm^2)');
ylabel('Capacitance (pF)');
title('Capacitance vs Plate Area');
grid on;

% Plot Capacitance vs Plate Separation
figure;
plot(sep_trials.("Plate Separation (mm)"), sep_trials.("Capacitance (pF)"), 's-', 'LineWidth', 2, 'Color', 'r');
xlabel('Plate Separation (mm)');
ylabel('Capacitance (pF)');
title('Capacitance vs Plate Separation');
grid on;
