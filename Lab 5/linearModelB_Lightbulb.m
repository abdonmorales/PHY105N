V=[5 6 7 8 9 10 11 12];
I = [0.698 0.758 0.822 0.873 0.928 0.978 1.025 1.071];

% Part I instrument specs
sigmaV = 0.01;     % V
sigmaI = 0.002;    % A

% --- Force column vectors ---
I = I(:); V = V(:);

% --- First, an unweighted quick fit to initialize m,b ---
X = [I, ones(size(I))];              % n×2
beta0 = X \ V;                       % [m0; b0]
m = beta0(1); b = beta0(2);

% --- Build effective vertical uncertainty: sigma_eff^2 = sigmaV^2 + (m*sigmaI)^2 ---
sigmaEff = sqrt( (sigmaV^2) * ones(size(V)) + (m*sigmaI).^2 );   % n×1

% --- Weighted least squares using lscov with weights vector w = 1/var ---
w = 1 ./ (sigmaEff.^2);              % n×1
beta  = lscov(X, V, w);              % [m; b], shape-safe
m = beta(1); b = beta(2);

% (Optional: one refinement pass)
sigmaEff = sqrt( (sigmaV^2) * ones(size(V)) + (m*sigmaI).^2 );
w = 1 ./ (sigmaEff.^2);
beta  = lscov(X, V, w);
m = beta(1); b = beta(2);

% --- Predictions, residuals, chi^2 ---
Vhat = X * beta;
res  = V - Vhat;
chi2 = sum( (res ./ sigmaEff).^2 );
N = numel(V); p = 2; dof = N - p;
redchi2 = chi2 / dof;
pval = 1 - chi2cdf(chi2, dof);

% --- Parameter uncertainties from covariance (shape-safe) ---
Cov = inv(X' * (w .* X));            % uses broadcasting; w is n×1
se  = sqrt(diag(Cov));
sm = se(1); sb = se(2);

fprintf('Model V = m I + b\n');
fprintf('m = %.6g ± %.2g ohms\n', m, sm);
fprintf('b = %.6g ± %.2g V\n', b, sb);
fprintf('chi2 = %.3f, dof = %d, red chi2 = %.3f, p = %.3f\n', chi2, dof, redchi2, pval);

% --- Plots ---
figure; errorbar(I, V, sigmaEff, 'o', 'LineStyle','none'); hold on;
plot(I, Vhat, '-', 'LineWidth', 1.6); grid on;
xlabel('I (A)'); ylabel('V (V)'); title('Lightbulb: V = m I + b');
legend('data (±σ_{eff})','fit','Location','best');

figure; plot(I, res, 'o-'); yline(0,'--'); grid on;
xlabel('I (A)'); ylabel('Residual (V - \hat{V})'); title('Residuals: V = m I + b');