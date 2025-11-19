%% ---------------- Single-slit: y vs (L/a) with error bars ----------------
% Model (lab): y = λ * (L/a) , intercept fixed at 0
% y in mm  -> λ returned in mm (convert to nm with *1e6)

% --- Data from your table ---
L_cm  = [74.5 70.9 45.5 51.1 76.1 56.2 74.2 49.5]';
a_cm  = [0.002 0.004 0.004 0.004 0.008 0.008 0.016 0.004]';
y_mm  = [1000 700 500 1600 30 300 200 500]';
x     = L_cm ./ a_cm;     % L/a (dimensionless)

% --- Uncertainties ---
sigma_a_cm = 0.001;       % ±0.01 mm = ±0.001 cm  (from your note)
% Propagate to x = L/a:  σx = |dx/da| σa = (L/a^2) * σa   (assuming σL negligible)
sigma_x = (L_cm ./ (a_cm.^2)) * sigma_a_cm;   % dimensionless

% y uncertainty (pick your ruler resolution; edit if you know it)
sigma_y_mm = 1.0 * ones(size(y_mm));  % e.g., ±1 mm

% --- Unweighted fit through the origin (what the lab expects) ---
lambda_mm = x \ y_mm;                    % slope in mm
lambda_nm = lambda_mm * 1e6;

% --- Optional: weighted fit through origin (uses σy) ---
W = diag(1 ./ (sigma_y_mm.^2));          % weights = 1/σ_y^2
lambda_w_mm = (x' * W * y_mm) / (x' * W * x);
lambda_w_nm = lambda_w_mm * 1e6;

% --- Compare to nominal laser with ±20 nm ---
lambda_nom_nm = 650;                     % set to your laser color if different
laser_tol_nm  = 20;
fprintf('Unweighted  λ = %.1f nm\n', lambda_nm);
fprintf('Weighted    λ = %.1f nm\n', lambda_w_nm);
fprintf('Laser nominal = %d ± %d nm\n\n', lambda_nom_nm, laser_tol_nm);

% --- Plot with error bars ---
figure; hold on;
% vertical (y) error bars:
errorbar(x, y_mm, sigma_y_mm, 'LineStyle','none', 'CapSize', 8);
% horizontal (x) error bars (drawn manually for wide MATLAB support)
for i = 1:numel(x)
    line([x(i)-sigma_x(i), x(i)+sigma_x(i)], [y_mm(i), y_mm(i)]);
end
scatter(x, y_mm, 36, 'filled');

% best-fit line (through origin)
xfit = linspace(min(x), max(x), 200)';
yfit = lambda_mm * xfit;
plot(xfit, yfit, 'LineWidth', 1.8);

xlabel('L/a (dimensionless)');
ylabel('y_{dark} (mm)');
title('y_{dark} vs (L/a): slope = \lambda (intercept = 0)');
legend('σ_y (vertical)','σ_x (horizontal)','Data','Fit','Location','best');

% Annotate result
txt = sprintf('\\lambda = %.4f mm  (%.0f nm)\\n(weighted: %.0f nm)', ...
              lambda_mm, lambda_nm, lambda_w_nm);
xpos = min(x) + 0.05*(max(x)-min(x));
ypos = min(y_mm) + 0.85*(max(y_mm)-min(y_mm));
text(xpos, ypos, txt, 'BackgroundColor', 'w', 'Margin', 4);
