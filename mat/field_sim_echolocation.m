%% Echolocation Simulation Across Velocities + Analysis
clear; clc; 
% close all;
savePlot = 0;
% =================== PARAMETERS ===================
fs = 192e3;                        % Sampling rate
bandwidth = [35e3, 65e3];          % Call bandwidth
kr = 5;                            % Responsivity coefficient
target_distance = 20;              % Initial distance (m)
initial_call_duration = 0.01;     % s
motile = true;                     % add motile noise
makeAudio = false;                 % audio off to save memory
Cr_max = 200;                      % buzz ceiling

% =================== STORAGE ===================
allRates = [];
allVels  = [];
allIPI   = [];
allDist  = [];
allSeq   = [];

% =================== RUN SIMULATIONS ===================
velocities = 1:0.05:20;   % range of initial velocities
allRows = [];            % accumulator for full long-format table
seqNum = 1;

for v0 = velocities
    fprintf('Running velocity %.2f m/s...\n', v0);

    % Run simulation
    result = simulateEcholocation(bandwidth, kr, target_distance, ...
                                  v0, initial_call_duration, ...
                                  motile, makeAudio);

    % --- Extract vectors ---
    Cr   = min(1 ./ result.delta_t(2:end), result.max_call_rate); % call rate (Hz)
    vr   = abs(result.delta_s ./ result.delta_t(2:end));          % instantaneous velocity
    d    = result.new_target_distance(1:end-1);                   % distance
    ipi  = result.delta_t(2:end);                                 % inter-pulse intervals
    dur  = result.new_call_duration(1:end-1);                     % call durations
    amp  = result.new_amplitude(1:end-1);                         % amplitudes

    % Validate equal lengths
    n = min([length(Cr), length(vr), length(d), length(ipi), length(dur), length(amp)]);
    Cr = Cr(1:n); vr = vr(1:n); d = d(1:n); ipi = ipi(1:n); dur = dur(1:n); amp = amp(1:n);

    % --- Repeat scalar values for each row ---
    row.SeqNum             = repmat(seqNum, n, 1);
    row.InitialVelocity    = repmat(v0, n, 1);
    row.Kr                 = repmat(result.Kr, n, 1);
    row.InitialTargetDist  = repmat(result.initial_target_distance, n, 1);
    row.Tb                 = repmat(result.Tb, n, 1);
    row.Tb_prime           = repmat(result.Tb_prime, n, 1);
    row.Tb_idx             = repmat(result.Tb_idx, n, 1);
    row.decrement_point    = repmat(result.decrement_point, n, 1);
    row.detection_point    = repmat(result.detection_point, n, 1);
    row.MaxCallRatePoint   = repmat(result.MaxCallRatePoint, n, 1);
    row.OverlapPoint       = repmat(result.OverlapPoint, n, 1);

    % --- Attach vectors ---
    row.Rate       = Cr(:);
    row.Velocity   = vr(:);
    row.Distance   = d(:);
    row.IPI        = ipi(:);
    row.CallDur    = dur(:);
    row.Amplitude  = amp(:);

    % --- Convert struct to table and append ---
    T = struct2table(row);
    allRows = [allRows; T]; %#ok<AGROW>

    seqNum = seqNum + 1;
end

% =================== SAVE TO CSV ===================
writetable(allRows, 'simulation_results.csv');
fprintf('Simulation results saved to simulation_results.csv\n');

%% =================== ANALYSIS PLOTS FOR SIMULATED DATA ===================

% Load simulation results
SimTable = readtable('simulation_results.csv');

% Unique sequences
seqs = unique(SimTable.SeqNum);

% Storage
allRates = [];
allVels  = [];
allIPI   = [];
allDist  = [];
allDur   = [];
allDur_ms = [];
allDurDiff  = [];
allDistDiff = [];

% --------- Process sequence by sequence ---------
for s = seqs'
    idx = SimTable.SeqNum == s;

    rates = SimTable.Rate(idx);
    vels  = SimTable.Velocity(idx);
    ipi   = SimTable.IPI(idx);
    dist  = SimTable.Distance(idx);
    dur   = SimTable.CallDur(idx);
    dur_ms = dur * 1000;

    % Differences within this sequence only
    durDiff  = [NaN; diff(dur)];
    distDiff = [NaN; diff(dist)];

    % Append
    allRates     = [allRates; rates];
    allVels      = [allVels; vels];
    allIPI       = [allIPI; ipi];
    allDist      = [allDist; dist];
    allDur       = [allDur; dur];
    allDur_ms    = [allDur_ms; dur_ms];
    allDurDiff   = [allDurDiff; durDiff];
    allDistDiff  = [allDistDiff; distDiff];
end

% ================== FILTER EXTREME VALUES ==================
validIdx = ...
    allRates > 0 & allRates < 200 & ...        % avoid wall at 200 Hz
    allDur_ms > 0.5 & allDur_ms < 10 & ...     % avoid min/max call duration walls
    allIPI > 1e-4 & allIPI < 0.5 & ...         % avoid capped IPIs
    allVels > 0 & allVels < 50;                % reasonable velocities

% Apply filter
allRates    = allRates(validIdx);
allVels     = allVels(validIdx);
allIPI      = allIPI(validIdx);
allDist     = allDist(validIdx);
allDur_ms   = allDur_ms(validIdx);
allDur      = allDur(validIdx);
allDurDiff  = allDurDiff(validIdx);
allDistDiff = allDistDiff(validIdx);

% ================== SUBSAMPLE FOR PLOTTING ==================
rng(1); % reproducible random subset
Nsample = 12000;
idx = randperm(length(allRates), min(Nsample, length(allRates)));

allRates    = allRates(idx);
allVels     = allVels(idx);
allIPI      = allIPI(idx);
allDist     = allDist(idx);
allDur_ms   = allDur_ms(idx);
allDurDiff  = allDurDiff(idx);
allDistDiff = allDistDiff(idx);
%% Figure with 6 panels
figure('Position', [100 100 500 500]);

% (a) Velocity vs Distance Covered
subplot(3,2,1);
x = abs(allDistDiff); y = allVels;
[r, p] = corr(y, x, 'Type','Spearman','Rows','complete');
scatter(x, y, 5, 'k','filled'); hold on;
coeffs = polyfit(x, y, 1);
xfit = linspace(min(x), max(x), 200);
plot(xfit, polyval(coeffs, xfit), 'r-','LineWidth',1.5);
xlabel('Distance Covered (m)','Interpreter','latex');
ylabel('Velocity (m/s)','Interpreter','latex');
title('$\mathbf{(a)}$', 'Interpreter','latex');
subtitle(sprintf('Velocity vs. Distance (r=%.3f, p=%.4f)', r, p), 'Interpreter','latex');
grid on; hold off;
axis square

% (b) Log Call Duration vs Velocity
subplot(3,2,2);
logDur = log10(allDur_ms./1000);
logVel = log10(allVels+eps);
[r, p] = corr(logDur, logVel, 'Type','Spearman','Rows','complete');
scatter(logVel, logDur, 5, 'k','filled'); hold on;
coeffs = polyfit(logVel, logDur, 1);
plot(sort(logVel), polyval(coeffs, sort(logVel)), 'r-','LineWidth',1.5);
xlabel('$log_{10}$ Velocity (m/s)','Interpreter','latex');
ylabel('$log_{10}$ Call Duration (ms)','Interpreter','latex');
title('$\mathbf{(b)}$', 'Interpreter','latex');
subtitle(sprintf('Log Call Duration vs Velocity (r=%.3f, p=%.4f)', r, p), 'Interpreter','latex');
xlim([-1 1.5])
grid on; axis square

% (c) Call Rate vs Call Duration
subplot(3,2,3);
[r, p] = corr(allRates, allDur_ms, 'Type','Spearman','Rows','complete');
scatter(allRates, allDur_ms, 5, 'k','filled');
xlabel('Call Rate (Hz)','Interpreter','latex');
ylabel('Call Duration (ms)','Interpreter','latex');
title('$\mathbf{(c)}$', 'Interpreter','latex');
subtitle(sprintf('Call Rate vs Call Duration (r=%.3f, p=%.4f)', r, p), 'Interpreter','latex');
grid on; axis square

% (d) Log Call Rate vs Log Call Duration
subplot(3,2,4);
logRate = log10(allRates+eps);
[r, p] = corr(logRate, logDur, 'Type','Spearman','Rows','complete');
scatter(logRate, logDur, 5, 'k','filled'); hold on;
coeffs = polyfit(logRate, logDur, 1);
plot(sort(logRate), polyval(coeffs, sort(logRate)), 'r-','LineWidth',1.5);
xlabel('$log_{10}$ Call Rate (Hz)','Interpreter','latex');
ylabel('$log_{10}$ Call Duration (ms)','Interpreter','latex');
title('$\mathbf{(d)}$', 'Interpreter','latex');
subtitle(sprintf('Log Call Rate vs Log Call Duration (r=%.3f, p=%.4f)', r, p), 'Interpreter','latex');
% ylim([-1 1])
xlim([0 2.5])
grid on; axis square

% (e) ΔCall Duration vs ΔDistance
subplot(3,2,5);
[r, p] = corr(allDurDiff, allDistDiff, 'Type','Spearman','Rows','complete');
scatter(allDistDiff, allDurDiff, 5, 'k','filled');
xlabel('Change in Distance Covered (m)','Interpreter','latex');
ylabel('Change in Call Duration (s)','Interpreter','latex');
title('$\mathbf{(e)}$', 'Interpreter','latex');
subtitle(sprintf('$\\Delta$Duration vs $\\Delta$Distance (r=%.3f, p=%.4f)', r, p), 'Interpreter','latex');
grid on; axis square
xlim([-2.5 2.5])
ylim([-0.015 0.01])

% (f) ΔCall Duration vs IPI
subplot(3,2,6);
[r, p] = corr(allDurDiff, allIPI, 'Type','Spearman','Rows','complete');
scatter(allDurDiff, allIPI*1000, 5, 'k','filled');
xlabel('Change in Call Duration (s)','Interpreter','latex');
ylabel('IPI (ms)','Interpreter','latex');
title('$\mathbf{(f)}$', 'Interpreter','latex');
subtitle(sprintf('$\\Delta$Duration vs IPI (r=%.3f, p=%.4f)', r, p), 'Interpreter','latex');
xlim([-0.015 0.015])
grid on; axis square

% Format layout
set(gcf,'Units','normalized','Position',[0.05 0.1 0.9 0.8]);
sgtitle('Simulated Echolocation Behaviour Profile','Interpreter','latex','FontSize',16);

%% Save figure
if savePlot
    saveFigure(gcf, '../biosonar_responsivity/fig', 'simulated_beh_profile')
end