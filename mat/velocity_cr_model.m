%% Multi-run Foraging Call Sequence with DBSCAN
% Compare Clean vs Noisy vs Field data using the *same* method
clear; close all;

% --- Parameters ---
c = 343;
kr = 5;
target_distance = 10; % m
max_call_rate = 200;  % calls/s
k_noise = 0.1;        % proportional tracking error

% --- Range of initial velocities ---
init_velocities = 1:0.5:50;   % m/s
rng(42);

% Storage
all_Cr_clean = []; all_vr_clean = []; all_d_clean = [];
all_Cr_noisy = []; all_vr_noisy = []; all_d_noisy = [];

for v0 = init_velocities
    % ----- Clean run -----
    new_target_distance(1) = target_distance;
    delta_t(1) = 0;
    call_number = 1; done = false;
    while ~done
        Ta = 2*new_target_distance(call_number)/c;
        delta_t(call_number+1) = (1+kr) * Ta;

        % step (clean)
        delta_s(call_number) = v0 * delta_t(call_number+1);

        new_target_distance(call_number+1) = ...
            new_target_distance(call_number) - delta_s(call_number);

        if new_target_distance(call_number) < 0.05, done = true; end
        call_number = call_number+1;
    end
    IPI = delta_t(2:end);
    Cr_i = min(1 ./ IPI, max_call_rate);
    vr_i = abs(delta_s ./ IPI);
    d_i  = new_target_distance(1:end-1);
    valid = isfinite(vr_i) & isfinite(Cr_i) & vr_i>0 & Cr_i>0;
    all_Cr_clean = [all_Cr_clean Cr_i(valid)];
    all_vr_clean = [all_vr_clean vr_i(valid)];
    all_d_clean  = [all_d_clean d_i(valid)];

    % ----- Noisy run -----
    new_target_distance(1) = target_distance;
    delta_t(1) = 0;
    call_number = 1; done = false;
    while ~done
        Ta = 2*new_target_distance(call_number)/c;
        delta_t(call_number+1) = kr * Ta;

        % step (with proportional noise)
        IPI_now = delta_t(call_number+1);
        eps_step = k_noise * (v0 * IPI_now) * randn();
        delta_s(call_number) = v0 * IPI_now + eps_step;

        new_target_distance(call_number+1) = ...
            new_target_distance(call_number) - delta_s(call_number);

        if new_target_distance(call_number) < 0.05, done = true; end
        call_number = call_number+1;
    end
    IPI = delta_t(2:end);
    Cr_i = min(1 ./ IPI, max_call_rate);
    vr_i = abs(delta_s ./ IPI);
    d_i  = new_target_distance(1:end-1);
    valid = isfinite(vr_i) & isfinite(Cr_i) & vr_i>0 & Cr_i>0;
    all_Cr_noisy = [all_Cr_noisy Cr_i(valid)];
    all_vr_noisy = [all_vr_noisy vr_i(valid)];
    all_d_noisy  = [all_d_noisy d_i(valid)];
end

% --- Apply same mask & DBSCAN + power-law fit ---
epsilon = 5; minpts = 6;
[Cr_clean, vr_clean, labels_clean, K_clean, b_clean, R2_clean] = ...
    filter_cluster_fit(all_Cr_clean, all_vr_clean, all_d_clean, max_call_rate, target_distance, epsilon, minpts, true);

[Cr_noisy, vr_noisy, labels_noisy, K_noisy, b_noisy, R2_noisy] = ...
    filter_cluster_fit(all_Cr_noisy, all_vr_noisy, all_d_noisy, max_call_rate, target_distance, epsilon, minpts, true);
% --- FIELD DATA ---
data = readtable('data/vof_processed_data.csv');
allRates = data.Rate; allVels = data.Velocity;
allDist  = 10 * ones(size(allRates));
[Cr_field, vr_field, labels_field, K_field, b_field, R2_field] = ...
    filter_cluster_fit(allRates, allVels, [], max_call_rate, target_distance, epsilon, minpts, false);

% --- Plotting ---
figure('Color','w', 'Position', [100 100 1000 400]); tiledlayout(1,3,'TileSpacing','compact');

% Clean
nexttile; hold on; grid on; box on;
gscatter(Cr_clean, vr_clean, labels_clean, lines(max(labels_clean)), '.', 10);
Cr_fit = logspace(log10(min(Cr_clean)), log10(max(Cr_clean)), 200);
plot(Cr_fit, K_clean * Cr_fit.^b_clean, 'm--','LineWidth',2);
xlabel('$C_r (Hz)$')
ylabel('$v~(m/s)$')
ylim([0 40])
title(sprintf('Clean Model\n$R^2=%.3f$, $v=%.2fC_r^{%.2f}$', R2_clean, K_clean, b_clean), 'Interpreter','latex');
axis square
legend off
formatLatex(gca)

% Noisy
nexttile; hold on; grid on; box on;
gscatter(Cr_noisy, vr_noisy, labels_noisy, lines(max(labels_noisy)), '.', 10);
Cr_fit = logspace(log10(min(Cr_noisy)), log10(max(Cr_noisy)), 200);
plot(Cr_fit, K_noisy * Cr_fit.^b_noisy, 'm--','LineWidth',2);
legend off
title(sprintf('Model with Noise\n$R^2=%.3f$, $v=%.2fC_r^{%.2f}$', R2_noisy, K_noisy, b_noisy), 'Interpreter','latex');
xlabel('$C_r (Hz)$')
ylabel('')
ylim([0 40])
axis square
formatLatex(gca)

% Field
nexttile; hold on; grid on; box on;
gscatter(Cr_field, vr_field, labels_field, lines(max(labels_field)), '.', 10);
Cr_fit = logspace(log10(min(Cr_field)), log10(max(Cr_field)), 200);
plot(Cr_fit, K_field * Cr_fit.^b_field, 'm--','LineWidth',2);
xlabel('$C_r (Hz)$')
ylabel('')
ylim([0 40])
legend off
title(sprintf('Field Data\n$R^2=%.3f$, $v=%.2fC_r^{%.2f}$', R2_field, K_field, b_field), 'Interpreter','latex');
axis square
formatLatex(gca)

sgtitle('Velocity Call-Rate Tradeoff ','Interpreter','latex', 'FontSize', 18);

fprintf('\n=== Power-law Fits ===\n');
fprintf('Clean: v_r = %.3f C_r^{%.3f}, R^2=%.4f\n', K_clean, b_clean, R2_clean);
fprintf('Noisy: v_r = %.3f C_r^{%.3f}, R^2=%.4f\n', K_noisy, b_noisy, R2_noisy);
fprintf('Field: v_r = %.3f C_r^{%.3f}, R^2=%.4f\n', K_field, b_field, R2_field);


%% Run bootstrap for each dataset
nboot = 1000;

[bCI_clean, KCI_clean] = bootstrap_powerlaw(Cr_clean, vr_clean, nboot);
[bCI_noisy, KCI_noisy] = bootstrap_powerlaw(Cr_noisy, vr_noisy, nboot);
[bCI_field, KCI_field] = bootstrap_powerlaw(Cr_field, vr_field, nboot);

% Print results
fprintf('\n=== Bootstrap 95%% Confidence Intervals ===\n');
fprintf('Clean Simulation: b = %.3f (95%% CI [%.3f, %.3f]), K = %.2f (95%% CI [%.2f, %.2f])\n', ...
    b_clean, bCI_clean(1), bCI_clean(2), K_clean, KCI_clean(1), KCI_clean(2));
fprintf('Noisy Simulation: b = %.3f (95%% CI [%.3f, %.3f]), K = %.2f (95%% CI [%.2f, %.2f])\n', ...
    b_noisy, bCI_noisy(1), bCI_noisy(2), K_noisy, KCI_noisy(1), KCI_noisy(2));
fprintf('Field Data:       b = %.3f (95%% CI [%.3f, %.3f]), K = %.2f (95%% CI [%.2f, %.2f])\n', ...
    b_field, bCI_field(1), bCI_field(2), K_field, KCI_field(1), KCI_field(2));

%% Save 
fig_path = '../biosonar_responsivity/fig';
% saveFigure(gcf, fig_path, 'vr_cr_model_data_fit.pdf')
%% --- Helper function ---
function [Cr_out, vr_out, labels_out, K, b, R2] = filter_cluster_fit(Cr_all, vr_all, d_all, maxCr, d0, epsDB, minptsDB, useDistance)
    % Filter points: always enforce valid call rates & velocities
    mask = (Cr_all < maxCr) & (Cr_all > 0) & (vr_all > 0);

    % Optionally apply distance filter (simulation only)
    if useDistance
        mask = mask & (d_all < 0.9 * d0);
    end

    Cr_use = Cr_all(mask); 
    vr_use = vr_all(mask);

    if isempty(Cr_use)
        Cr_out=[]; vr_out=[]; labels_out=[]; K=NaN; b=NaN; R2=NaN; labels_out=[];
        return;
    end

    % DBSCAN
    X = [Cr_use(:), vr_use(:)];
    labels = dbscan(X, epsDB, minptsDB);
    keep   = labels > 0;

    Cr_out = Cr_use(keep); 
    vr_out = vr_use(keep);
    labels_out = labels(keep);

    if numel(Cr_out) < 3
        K=NaN; b=NaN; R2=NaN; return;
    end

    % Power-law fit
    lx = log(Cr_out(:)); ly = log(vr_out(:));
    p  = polyfit(lx, ly, 1);
    b  = p(1); K = exp(p(2));
    lyhat = polyval(p, lx);
    R2 = 1 - sum((ly-lyhat).^2)/sum((ly-mean(ly)).^2);
end

% --- Bootstrapping confidence intervals ---
function [b_CI, K_CI] = bootstrap_powerlaw(Cr, vr, nboot)
    lx = log(Cr(:)); ly = log(vr(:));
    n = numel(lx);
    boot_b = zeros(nboot,1);
    boot_K = zeros(nboot,1);
    for i = 1:nboot
        idx = randi(n, n, 1);   % resample indices
        p = polyfit(lx(idx), ly(idx), 1);
        boot_b(i) = p(1);
        boot_K(i) = exp(p(2));
    end
    b_CI = prctile(boot_b, [2.5 97.5]);
    K_CI = prctile(boot_K, [2.5 97.5]);
end