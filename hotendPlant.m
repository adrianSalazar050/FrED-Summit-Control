%% import Data

data = importdata('heaterData.txt', ',');
t = data(:,1)/1000;
T = data(:,3);

T0 = mean(T(1:20));
y  = T - T0;
y_ss = mean(y(end-200:end));

%% Initial estimates
L_est   = t(find(y > 0.01*y_ss,1));
tau_est = t(find(y >= 0.632*y_ss,1)) - L_est;

% Model tf
fopdt_tf = @(p,t) p(1) * step( tf(1,[p(2) 1],'InputDelay',p(3)), t );

% Fit
p = lsqcurvefit(fopdt_tf, [y_ss tau_est L_est], t, y, [0 10 0], [y_ss*1.5 2000 60]);
K = p(1); tau = p(2); L = p(3);

% Continuous transfer function
G = tf(K, [tau 1], 'InputDelay', L);

%% Discretization (Ts = 0.12 s = 120 ms)
Ts = 0.12;
Gd = c2d(G, Ts, 'zoh');

%% Simulated response
y_model = step(G, t);

%% R^2
ss_res = sum((y - y_model).^2);
ss_tot = sum((y - mean(y)).^2);
R2 = 1 - ss_res/ss_tot;

fprintf('R^2 = %.6f\n', R2);

% Plot
figure
plot(t/60, T, 'b'); hold on
plot(t/60, y_model + T0, '--r','LineWidth',2)
grid on
xlabel('Time (min)')
ylabel('Temperature (°C)')
legend('Measured', sprintf('FOPDT (R^2 = %.4f)', R2))
title('Simplified FOPDT Fit')