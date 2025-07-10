% Plot resonsivity curves and inflection poin demo
 initial_call_duration = 0.005;
 motile = 0;
 makeAudio = 0;
bandwidth = 65e3:-100:35e3; 
target_distance = 10;            % Initial target distance (meters)
initial_velocity = 5;
Cr_max = 200;
result = simulateEcholocation(bandwidth, initial_velocity, target_distance, 3, initial_call_duration, motile, makeAudio);

R = 1./diff(result.delta_t);
Rmax = abs(abs(R)-Cr_max);

% === Data Setup ===
delta_t = result.delta_t;
R = 1 ./ diff(delta_t);
Cr_max = 200;
[~, Tb_idx] = min(abs(abs(R) - Cr_max));
ddt = diff(delta_t);
Tb_prime = ddt(Tb_idx);
Rmax = abs(abs(R) - Cr_max);

% === Plot Side-by-Side ===
figure('Color','w', 'Position', [100, 100, 1200, 500]);

% --- Subplot 1: Responsivity Analysis ---
subplot(1,2,1);
plot(abs(R), 'b--', 'LineWidth', 2); hold on;
plot(Rmax, 'm-', 'LineWidth', 2);
plot(Tb_idx, Rmax(Tb_idx), 'k+', 'MarkerSize', 20);
ylabel("Responsivity, $\mathcal{R}$", 'Interpreter', 'latex');
xlabel("Call Number");
xt = xticks;
xlim([0 max(xt)]);
axis square
formatLatex(gca);

% Add Responsivity equation
annotation('textbox', [0.22, 0.75, 0.25, 0.1], ...
    'String', '$$\mathcal{R}_n = \left|\frac{1}{\Delta t_{n+1} - \Delta t_n}\right|$$', ...
    'Interpreter', 'latex', 'FontSize', 14, 'Color', 'b', 'EdgeColor', 'none');

% Add argmin equation
annotation('textbox', [0.22, 0.65, 0.25, 0.1], ...
    'String', '$$n^* = \arg\min_n |\mathcal{R}_n - C_{r,\mathrm{max}}|$$', ...
    'Interpreter', 'latex', 'FontSize', 14, 'Color', 'm', 'EdgeColor', 'none');

% --- Subplot 2: Tb* from IPI Difference ---
subplot(1,2,2); hold on;
plot(1:length(delta_t), delta_t*1000, '-o', 'MarkerFaceColor', 'b');
xlabel('Call Number');
ylabel('$\Delta t$ (IPI), ms', 'Interpreter', 'latex');
title('$T_{b^*}$ from $\mathcal{R}$ \& IPIs', 'Interpreter','latex');
xlim([2 max(xt)]);
axis square
% Mark the Tb* computation
plot([Tb_idx Tb_idx+1], delta_t([Tb_idx Tb_idx+1])*1000, 'go', 'MarkerFaceColor', 'm');
line([Tb_idx+0.5 Tb_idx+0.5], delta_t([Tb_idx Tb_idx+1])*1000, ...
    'Color', 'g', 'LineStyle', '--', 'LineWidth', 2);
text(Tb_idx + 1, mean(delta_t([Tb_idx Tb_idx+1])*1000), ...
     sprintf('$T_{b^*} = %.2f$ ms', abs(1000*Tb_prime)), ...
     'Interpreter', 'latex', 'FontSize', 14, 'Color', 'm');

formatLatex(gca);
%%
img_path = '/Users/ravi/Documents/projects/thesis/papers/biosonar_responsivity/fig/R_curve';
exportgraphics(gcf, [img_path '.pdf'], 'Resolution', 300, 'Append', false)