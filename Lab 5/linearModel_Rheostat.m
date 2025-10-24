V=[5 6 7 8 9 10 11 12];
I = [0.032 0.038 0.044 0.051 0.057 0.064 0.069 0.075];
sigmaV = 0.01;
% --- Make column vectors and weights ---
I = I(:); V = V(:);
sigmaV = sigmaV * ones(size(V));      % n×1
w = 1./(sigmaV.^2);                   % n×1

% ===== MODEL 1: V = m I (through origin) =====
m = (I'*(w.*V)) / (I'*(w.*I));        % weighted slope
Vhat = m*I; res = V - Vhat;

chi2 = sum((res./sigmaV).^2);
n = numel(V); p = 1; dof = n - p;
redchi2 = chi2/dof;  pval = 1 - chi2cdf(chi2,dof);

fprintf('V = mI:  m = %.6g ohms\n', m);
fprintf('chi2 = %.3f, dof = %d, red chi2 = %.3f, p = %.3f\n', chi2, dof, redchi2, pval);

% Plot
figure; errorbar(I,V,sigmaV,'o','LineStyle','none'); hold on; plot(I,Vhat,'-');
xlabel('I (A)'); ylabel('V (V)'); title('V = mI'); grid on; legend('data (±0.01 V)','fit');