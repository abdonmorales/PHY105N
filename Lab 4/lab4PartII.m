%% Part II — raw data from your sheet (no CSV)

% Wavelength (given on sheet)
lambda  = 530e-9;     % m
dlambda = 20e-9;      % m

% -------- HOLES (trials 1–8) --------
L_holes_cm = [79.9 62.4 82.1 57.3 44.0 83.9 31.6 38.2];
a_holes_mm = 0.6 * ones(1,8);
y_holes_m  = [1 0.8 0.9 0.6 0.5 0.9 0.4 0.4];  % already in meters per your header
scale_holes = ones(1,8);                       % change if you meant "ydark scale"

% -------- CIRCULAR (trials 9–17) --------
L_circ_cm = [83.4 65.9 43.2 73.9 78.6 70.9 76.5 82.8 55.1];
a_circ_mm = [0.2 0.2 0.2 0.2 0.2 0.2 0.4 0.4 0.4];
y_circ_m  = [0.10 0.15 1.00 0.12 0.009 0.045 0.00021 0.00020 0.00004];
scale_circ = ones(1,9);                        % change if needed

%% ---- Helper: zero-intercept fit + 1σ on slope ----
zeroInterceptFit = @(xk,yk) deal( ...
    (xk'*yk)/(xk'*xk), ...                               % slope m
    sqrt( sum((yk - ((xk'*yk)/(xk'*xk)).*xk).^2) / max(numel(xk)-1,1) / sum(xk.^2) ) ...
);

%% ---- Pack patterns into a struct array to loop cleanly ----
patterns(1).name = "holes";
patterns(1).L_m  = L_holes_cm/100;          % cm → m
patterns(1).a_m  = a_holes_mm/1000;         % mm → m
patterns(1).y_m  = y_holes_m .* scale_holes;

patterns(2).name = "circular";
patterns(2).L_m  = L_circ_cm/100;
patterns(2).a_m  = a_circ_mm/1000;
patterns(2).y_m  = y_circ_m .* scale_circ;

%% ---- Loop over patterns: compute x=L/a, fit, and plot ----
for k = 1:numel(patterns)
    p = patterns(k);
    x = p.L_m ./ p.a_m;          % dimensionless
    y = p.y_m;                   % meters

    % Fit with intercept constrained to 0
    [m, dm] = zeroInterceptFit(x(:), y(:));  % slope m (m), 1σ dm (m)
    g = m / lambda;

    % Plot
    figure('Name',char(p.name)); hold on; grid on;
    scatter(x, y, 50, 'filled', 'DisplayName','Data');
    xlabel('L/a (dimensionless)'); ylabel('\Delta y (m)');
    title(sprintf('Part II: %s (\\Delta y vs L/a, intercept = 0)', p.name));

    xx = linspace(min(x), max(x), 200);
    plot(xx, m*xx, 'LineWidth', 1.8, 'DisplayName', ...
         sprintf('Fit: m = %.3e \\pm %.1e m', m, dm));

    % Reference model lines: lambda and lambda ± 20 nm
    plot(xx,  lambda*xx,  '--', 'LineWidth', 1.0, 'DisplayName','Model: \lambda = 530 nm');
    plot(xx, (lambda+dlambda)*xx, ':',  'LineWidth', 1.0, 'DisplayName','\lambda + 20 nm');
    plot(xx, (lambda-dlambda)*xx, ':',  'LineWidth', 1.0, 'DisplayName','\lambda - 20 nm');

    legend('Location','best');

    % Console summary
    fprintf('\nPattern: %s\n', p.name);
    fprintf('  slope m = %.3e m  (± %.1e)\n', m, dm);
    fprintf('  g = m/λ = %.3f   (λ = 5.30e-7 m)\n', g);
end