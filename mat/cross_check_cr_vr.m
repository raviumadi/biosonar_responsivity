%% Cross-check responsivity identities (discrete product form)
clear; 
cls;
savePlot = 0;
% ------------ Run simulation ------------
fs      = 192e3;
bw      = [45e3 90e3];
kr      = 5;              % responsivity coefficient
d0      = 10;             % m
v_set   = 5.0;            % m/s (toward target)
tcall0  = 0.005;          % s
motile  = false;
makeAud = false;

res = simulateEcholocation(bw, kr, d0, v_set, tcall0, motile, makeAud);

% ------------ Extract arrays ------------
c    = 343;                            % m/s
dt   = res.delta_t(:);                 % length N+1 with dt(1)=0
dt   = dt(2:end);                      % length N (per-step IPI)
Cr   = 1 ./ dt;                        % per-step call rate (Hz), tied to distance BEFORE each step

d    = res.new_target_distance(:);     % length N+1: d(1)=initial, d(n+1)=after step n
d_before = d(1:end-1);                 % length N: distance used to compute dt (and Cr)
d_after  = d(2:end);                   % length N: distance *after* the displacement

% Optional time stamps at call onsets (seconds)
t_calls = res.call_points(:) ./ fs;    % length N+1
t_step  = t_calls(2:end);              % time at end of each step (aligns with dt/d_after)

% ------------ Identities / theory ------------
% 1) Distance-from-call-rate identity (must match d BEFORE the step)
d_fromCr = c ./ ( 2*(1+kr) .* Cr );    % length N, compare to d_before

% 2) True relative speed from simulator (per step)
v_true = res.delta_s(:) ./ dt;         % length N, equals v_set when motile=false

% 3) Discrete product form for velocity:
% ΔCr ≈ (∂Cr/∂d)*Δd  ⇒  v_r ≈ c/[2(1+kr)] * (ΔCr/Δt) / Cr^2
dCr    = diff(Cr);                     % length N-1
dCr_dt = dCr ./ dt(1:end-1);             % length N-1
Cr_mid = 0.5 * (Cr(1:end-1) + Cr(2:end));
t_mid  = 0.5 * (t_step(1:end-1) + t_step(2:end));

v_theory = (c ./ (2*(1+kr))) .* (dCr_dt ./ (Cr_mid.^2));  % length N-1

% Midpoint version of "true" speed for apples-to-apples comparison
v_true_mid = 0.5 * (v_true(1:end-1) + v_true(2:end));     % length N-1

% ------------ Plots ------------
figure('Color','w', 'Position', [100 100 800 800]);

% (1) Velocity consistency (discrete-corrected)
subplot(2,2,1);
plot(t_mid, v_true_mid, 'k-', 'LineWidth', 1.5); hold on;
plot(t_mid, v_theory,   'r--','LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Relative velocity $v_r$ (m/s)');
legend('From simulator step $(\Delta s/\Delta t)$','From discrete call-rate product',...
    'Location','best', 'Interpreter', 'latex');
title('$\textbf{I.}$ Relative velocity consistency');
ylim([2 8])
formatLatex(gca);

% (2) Distance identity: BEFORE-step distance vs c/[2(1+k_r) C_r]
subplot(2,2,2);
plot(t_calls(1:end-1), d_before, 'k-', 'LineWidth', 1.5); hold on;
plot(t_calls(1:end-1), d_fromCr, 'b--','LineWidth', 2);
xlabel('Time (s)'); ylabel('Distance to target d (m)');
legend('Simulated d (before step)','$d = c/(2(1+k_r)C_r)$','Location','best',...
    'Interpreter', 'latex');
title('$\textbf{II.}$ Distance identity');
formatLatex(gca);

% (3) C_r – d (should lie on a curve)
subplot(2,2,3);
plot(d_before, Cr, 'k.', 'MarkerSize', 20); hold on;
plot(d_fromCr, Cr, 'r.', 'MarkerSize', 10);
xlabel('Distance d (m)'); ylabel('Call rate $C_r$ (Hz)');
title('$\textbf{III.}$ $C_r(d) = c / [2(1+k_r) d]$');
legend('Simulated pairs ($d,C_r)$','Identity using $C_r$','Location','best',...
    'Interpreter', 'latex');
ylim([0 220])
formatLatex(gca);

% (4) v_r vs C_r
subplot(2,2,4);
plot(Cr_mid, v_true_mid, 'k.', 'MarkerSize', 8); hold on;
plot(Cr_mid, v_theory,   'ro','MarkerSize', 10);
xlabel('Call rate $C_r$ (Hz)'); ylabel('Relative velocity $v_r$ (m/s)');
title('$\textbf{IV.}$ $v_r$ vs $C_r$');
legend('From simulator','From discrete theory','Location','best',...
    'Interpreter', 'latex');
ylim([2 8])
xlim([0 220])
formatLatex(gca);

sgtitle('Responsivity Framework: Validating Derived Parameters', 'Fontsize', 18,'Interpreter', 'latex' );

% ------------ Sanity metrics ------------
fprintf('Mean |v_true_mid - v_theory|: %.4f m/s\n', mean(abs(v_true_mid - v_theory)));
fprintf('Median |d_before - d_fromCr|: %.4f m\n', median(abs(d_before - d_fromCr)));


%% Save figure
if savePlot
    saveFigure(gcf, '/Users/ravi/Documents/projects/thesis/papers/biosonar_responsivity/fig', 'equation_validation')
end