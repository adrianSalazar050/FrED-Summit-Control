data = importdata('heaterData.txt', ',');
t    = data(:,1) / 1000;
T    = data(:,3);
Ts   = t(2) - t(1);

T0   = mean(T(1:20));
y    = T - T0;
y_ss = mean(y(end-200:end));

L_est   = t(find(y > 0.01 * y_ss, 1));
tau_est = t(find(y >= 0.632 * y_ss, 1)) - L_est;

fopdt = @(p,t) p(1) .* (t > p(3)) .* (1 - exp(-max(t-p(3),0)/p(2)));
opts = optimoptions('lsqcurvefit','Algorithm','trust-region-reflective','MaxIterations',2000,'Display','off');
p_fo = lsqcurvefit(fopdt, [y_ss, tau_est, L_est], t, y, [0, 10, 0], [y_ss*1.5, 2000, 60], opts);
K_fo = p_fo(1);  tau_fo = p_fo(2);  L_fo = p_fo(3);
y_fo = fopdt(p_fo, t);

ss_res = sum((y - y_fo).^2);
r2     = 1 - ss_res / sum((y - mean(y)).^2);
fprintf('FOPDT:  K=%.4f°C  τ=%.2fs  L=%.2fs  R²=%.6f\n', K_fo, tau_fo, L_fo, r2);

G_fo = tf(K_fo, [tau_fo, 1], 'InputDelay', L_fo);
fprintf('\nFOPDT Transfer Function G(s):\n'); disp(G_fo)

figure('Name','Hotend System ID','Position',[100 100 900 450]);
plot(t/60, T, 'Color',[0.12 0.47 0.71], 'LineWidth',1.2); hold on
plot(t/60, y_fo + T0, '--r', 'LineWidth',2)
xlabel('Time (min)'); ylabel('Temperature (°C)');
title('Hotend Step Response — FOPDT Fit','FontWeight','bold')
legend('Measured', sprintf('FOPDT  K=%.1f°C  τ=%.1fs  L=%.1fs  R²=%.4f', ...
       K_fo, tau_fo, L_fo, r2), 'Location','southeast')
grid on; xlim([0 t(end)/60])