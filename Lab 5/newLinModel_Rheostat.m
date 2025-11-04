V=[5;6;7;8;9;10;11;12];
I = [0.032;0.038;0.044;0.051;0.057;0.064;0.069;0.075];

m = 158.312;                  % from your fit
sigmaV = 0.01;                % V (Part I)
sigmaI = 0.002;               % A (Part I)

% effective vertical uncertainty used for chi^2
sigmaEff = sqrt( (sigmaV^2) + (m*sigmaI)^2 ) * ones(size(V));

Vhat = m*I;
res  = V - Vhat;

% Main plot
figure; 
errorbar(I, V, sigmaEff, 'o', 'LineStyle','none', 'MarkerSize', 6); hold on;
plot(I, Vhat, '-', 'LineWidth', 1.6);
xlabel('Current I (A)'); ylabel('Voltage V (V)');
title(sprintf('Rheostat: V = m I  (m = %.3f \\Omega)', m));
legend('data (\pm\sigma_{eff})','fit','Location','best'); grid on;