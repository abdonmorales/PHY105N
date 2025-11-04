xr = [2.928523524;2.809402695;2.564949357;2.501435952;2.312535424;2.251291799;2.208274414;2.104134154;2.041220329;1.987874348];
yF = [4.585987367;5.279134547;5.684599655;5.972281728;6.195425279;6.377746836;6.531897516;6.665428908;6.783211944;6.88857246];
sigma_y=[];

xr = xr(:); yF = yF(:);
n = numel(xr);
if numel(yF) ~= n, error('xr and yF must have same length'); end

%% === FIT: unweighted or weighted (if sigma_y provided) ===
if isempty(sigma_y)
    % Unweighted OLS
    p = polyfit(xr, yF, 1);            % yF = alpha*xr + B
    alpha = p(1);  B = p(2);
    yhat = polyval(p, xr);
    res  = yF - yhat;
    s2   = sum(res.^2) / (n - 2);      % residual variance
    Sxx  = sum( (xr - mean(xr)).^2 );
    alpha_stderr = sqrt(s2 / Sxx);
    % --- Chi-square style outputs for unweighted OLS (no given sigmas) ---
    RSS = sum(res.^2);          % residual sum of squares
    dof = n - 2;                % two fit params: alpha, B
    chi2 = RSS / s2;            % = sum( (res.^2) / s2 )
    chi2_red = chi2 / dof;      % reduced chi^2
    fprintf(['\nalpha = %.4f ± %.4f (SE)\nB = %.4f\nR^2 = %.4f\n' ...
         'RSS = %.6f\ns^2 = %.6f\ndof = %d\nchi^2 = %.6f\nchi^2_red = %.6f\n'], ...
         alpha, alpha_stderr, B, R2, RSS, s2, dof, chi2, chi2_red);

else
    % Weighted least squares with weights = 1/sigma_y^2
    sigma_y = sigma_y(:);
    if numel(sigma_y) ~= n, error('sigma_y must match length of xr/yF'); end
    w = 1 ./ (sigma_y.^2);
    X = [xr ones(n,1)];
    % (X' W X) beta = X' W y
    W = diag(w);
    beta = (X' * W * X) \ (X' * W * yF);
    alpha = beta(1); B = beta(2);
    yhat = X*beta; res = yF - yhat;
    % Weighted residual variance with dof = n-2
    s2 = sum(w .* (res.^2)) / (n - 2);
    % Covariance matrix of beta
    CovB = s2 * inv(X' * W * X);
    alpha_stderr = sqrt(CovB(1,1));
end

R2 = 1 - sum((yF - yhat).^2) / sum((yF - mean(yF)).^2);
%%fprintf('\nalpha = %.4f ± %.4f (SE)\nB = %.4f\nR^2 = %.4f\n', alpha, alpha_stderr, B, R2);

%% === PLOTS ===
% 1) Linearized plot: y = alpha x + B
figure(1); clf;
plot(xr, yF, 'o', 'MarkerSize', 7, 'LineWidth', 1.5); hold on;
plot(xr, yhat, '-', 'LineWidth', 2);
grid on;
xlabel('ln r'); ylabel('ln F');
title('Linearized Fit: ln F = \alpha ln r + B');
legend('Data', sprintf('\\alpha = %.3f \\pm %.3f,  B = %.3f', alpha, alpha_stderr, B), 'Location', 'best');

% 2) Reconstructed log–log plot by exponentiating
r  = exp(xr);
F  = exp(yF);
rfit = linspace(min(r), max(r), 300);
Ffit = exp(B) .* (rfit.^alpha);

figure(2); clf;
loglog(r, F, 'o', 'MarkerSize', 7, 'LineWidth', 1.5); hold on;
loglog(rfit, Ffit, '-', 'LineWidth', 2);
grid on;
xlabel('Distance r (m)'); ylabel('Force F (N)');
title('Magnetic Force vs Distance (Power-Law Fit)');
legend('Data', sprintf('Fit: F = K r^{\\alpha},  \\alpha = %.3f \\pm %.3f', alpha, alpha_stderr), 'Location', 'best');
