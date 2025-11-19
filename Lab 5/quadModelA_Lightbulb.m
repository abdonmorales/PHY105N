V=[5;6;7;8;9;10;11;12];
I = [0.698;0.758;0.822;0.873;0.928;0.978;1.025;1.071];
sigmaV = 0.01;  sigmaI = 0.002;

% Design matrix for V = a I^2 + m I (+ b optional)
X = [I.^2, I];                % no intercept; use [I.^2, I, ones(size(I))] if you want +b

% 1) Unweighted init
beta = X \ V;                 % [a; m]  (or [a;m;b] if you added ones)
for k = 1:5                   % a few Gauss-Newton style refinements
    a = beta(1); m = beta(2);
    dVdI = 2*a*I + m;
    sigmaEff = sqrt( (sigmaV^2) + (dVdI.*sigmaI).^2 );
    w = 1./(sigmaEff.^2);     % n×1 weights
    % weighted least squares
    beta = (X'*(w.*X)) \ (X'*(w.*V));
end

a = beta(1); m = beta(2);
Vhat = X*beta;
res  = V - Vhat;

chi2 = sum((res./sigmaEff).^2);
N = numel(V); p = size(X,2); dof = N - p;
redchi2 = chi2/dof;  pval = 1 - chi2cdf(chi2, dof);

% parameter uncertainties
Cov = inv(X'*(w.*X));  se = sqrt(diag(Cov));

fprintf('Model V = a I^2 + m I\n');
fprintf('a = %.6g ± %.2g V/A^2\n', a, se(1));
fprintf('m = %.6g ± %.2g ohms\n',    m, se(2));
fprintf('chi2=%.3f, dof=%d, red chi2=%.3f, p=%.3f\n', chi2, dof, redchi2, pval);

% Plot + residuals
figure; errorbar(I,V,sigmaEff,'o','LineStyle','none'); hold on;
plot(I,Vhat,'-','LineWidth',1.5); grid on;
xlabel('I (A)'); ylabel('V (V)'); title('Lightbulb: V = a I^2 + m I');
legend('data (±σ_{eff})','fit','Location','best');