%% Explore the shortest behavioural interval T_b at d_stop
% This script focuses only on the minimum behavioural interval implied by
% the responsivity framework at the stop/interception distance d_stop.
%
% Step 1: For each k_r, compute the theoretical maximum call rate at
% d_stop with v_r = 0. This gives the corresponding minimum IPI and its
% decomposition into T_a and T_b just before interception under the
% stationary-anchor reference.
%
% Step 2: For a range of radial closing velocities, compute the shortest
% T_b at d_stop and show how it varies with k_r and v_r.
%
% Step 3: Define a near-terminal spatial interval as the portion of the
% approach between f*C_r,max and C_r,max, where C_r,max is referenced to
% v_r = 0 at d_stop. We then compute the time spent in that final spatial
% interval and the corresponding number of calls. Under this formulation,
% time is controlled mainly by v_r, whereas call density within that same
% interval is controlled by k_r.
%
% Step 4: Illustrate how different late-approach deceleration profiles
% change the time and call budget inside that same final spatial interval.
% This makes explicit how slowing before capture can lengthen the buzz
% without changing the underlying responsivity coefficient.

clear; clc;

thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, 'fcn'));
outDir = fullfile(thisDir, 'results', 'exploratory_figures');
savePlots = false;
runConfig = responsivityRunConfig();
if runConfig.OverrideOutputSwitches
    savePlots = runConfig.SaveFigures;
end
if runConfig.CloseFiguresBeforeRun
    close all;
end

if savePlots && ~exist(outDir, 'dir')
    mkdir(outDir);
end

%% Parameters
c = 343;                       % speed of sound, m/s
dStop_m = 0.15;                % stop/interception distance
krValues = [3 5 7 9];
vrRange_m_s = 0.1:0.1:3;
finalFraction = 0.75;
colorsKr = parula(numel(krValues));

%% Step 1: stationary-anchor reference at v_r = 0
Tref = table();
for i = 1:numel(krValues)
    kr = krValues(i);

    TaMin_s = 2 * dStop_m / c;
    IpiMin_s = (1 + kr) * TaMin_s;
    TbMin_s = kr * TaMin_s;
    CrMax_Hz = 1 / IpiMin_s;

    newRow = table(kr, CrMax_Hz, IpiMin_s, TaMin_s, TbMin_s, ...
        'VariableNames', {'kr','CrMax_vr0_Hz','IpiMin_vr0_s','TaMin_vr0_s','TbMin_vr0_s'});
    Tref = [Tref; newRow]; %#ok<AGROW>
end

disp('=== Stationary-anchor reference at d_stop ===');
disp(Tref);

%% Step 2: shortest T_b at d_stop across radial closing velocities
Tall = table();
for i = 1:numel(krValues)
    kr = krValues(i);

    for j = 1:numel(vrRange_m_s)
        vr = vrRange_m_s(j);

        TaMin_s = 2 * dStop_m / (c + vr);
        IpiMin_s = (1 + kr) * TaMin_s;
        TbMin_s = kr * TaMin_s;
        CrMax_Hz = 1 / IpiMin_s;

        newRow = table(kr, vr, CrMax_Hz, IpiMin_s, TaMin_s, TbMin_s, ...
            'VariableNames', {'kr','vr_m_s','CrMax_Hz','IpiMin_s','TaMin_s','TbMin_s'});
        Tall = [Tall; newRow]; %#ok<AGROW>
    end
end

fprintf('d_stop = %.3f m\n', dStop_m);
fprintf('Across v_r = %.1f--%.1f m/s, shortest T_b ranges from %.3f to %.3f ms\n', ...
    min(vrRange_m_s), max(vrRange_m_s), ...
    1000 * min(Tall.TbMin_s), 1000 * max(Tall.TbMin_s));

%% Step 3: time and call budget in the final spatial interval
Tterminal = table();
for i = 1:numel(krValues)
    kr = krValues(i);
    CrMaxRef_Hz = Tref.CrMax_vr0_Hz(Tref.kr == kr);

    for j = 1:numel(vrRange_m_s)
        vr = vrRange_m_s(j);

        dFraction_m = (c + vr) / (2 * (1 + kr) * finalFraction * CrMaxRef_Hz);
        finalSpatialInterval_m = max(dFraction_m - dStop_m, 0);
        finalIntervalDuration_s = finalSpatialInterval_m / vr;

        if dFraction_m > dStop_m
            nCallsFinalInterval = ((c + vr) / (2 * (1 + kr) * vr)) * log(dFraction_m / dStop_m);
        else
            nCallsFinalInterval = 0;
        end

        newRow = table(kr, vr, CrMaxRef_Hz, dFraction_m, finalSpatialInterval_m, ...
            finalIntervalDuration_s, nCallsFinalInterval, ...
            'VariableNames', {'kr','vr_m_s','CrMaxRef_Hz','dFraction_m', ...
            'FinalSpatialInterval_m','FinalIntervalDuration_s','NCallsFinalInterval'});
        Tterminal = [Tterminal; newRow]; %#ok<AGROW>
    end
end

fprintf('Using the final interval from %.0f%% C_r,max to d_stop:\n', 100 * finalFraction);
fprintf('Time spent in that interval ranges from %.3f to %.3f s\n', ...
    min(Tterminal.FinalIntervalDuration_s), max(Tterminal.FinalIntervalDuration_s));
fprintf('Call count in that interval ranges from %.2f to %.2f calls\n', ...
    min(Tterminal.NCallsFinalInterval), max(Tterminal.NCallsFinalInterval));

%% Step 4: example late-approach deceleration profiles in distance space
krExample = 5;
d0Scenario_m = 1.2;
CrMaxRefExample_Hz = Tref.CrMax_vr0_Hz(Tref.kr == krExample);
dFractionRef_m = dStop_m / finalFraction;
dPath_m = linspace(d0Scenario_m, dStop_m, 500);

profileNames = {'No deceleration', 'Mild deceleration', 'Strong deceleration', 'Very strong deceleration'};
profileColors = [
    0.20 0.20 0.20
    0.85 0.45 0.10
    0.10 0.45 0.75
    0.55 0.20 0.70
];
profileStyles = {'--', '-', '-', '-'};
vStart_m_s = 2.5;
vStopValues_m_s = [2.5, 1.2, 0.5, 0.1];

Tprofiles = table();
for i = 1:numel(vStopValues_m_s)
    vStop = vStopValues_m_s(i);
    frac = (d0Scenario_m - dPath_m) / (d0Scenario_m - dStop_m);
    frac = min(max(frac, 0), 1);
    vPath_m_s = vStart_m_s + (vStop - vStart_m_s) .* frac;

    finalRows = dPath_m <= dFractionRef_m & dPath_m >= dStop_m;
    dFinal_m = dPath_m(finalRows);
    vFinal_m_s = vPath_m_s(finalRows);

    % Integrate over decreasing distance; flip to ascending distance for trapz.
    dAsc_m = fliplr(dFinal_m);
    vAsc_m_s = fliplr(vFinal_m_s);
    CrAsc_Hz = (c + vAsc_m_s) ./ (2 * (1 + krExample) .* dAsc_m);

    finalDuration_s = trapz(dAsc_m, 1 ./ vAsc_m_s);
    finalCalls = trapz(dAsc_m, CrAsc_Hz ./ vAsc_m_s);
    if abs(vStart_m_s - vStop) < eps
        decel_m_s2 = 0;
    else
        decel_m_s2 = (vStart_m_s^2 - vStop^2) / (2 * (d0Scenario_m - dStop_m));
    end

    newRow = table(string(profileNames{i}), vStart_m_s, vStop, finalDuration_s, ...
        finalCalls, mean(vAsc_m_s, 'omitnan'), ...
        'VariableNames', {'Profile', 'VelocityAtStart_m_s', 'VelocityAtStop_m_s', ...
        'FinalIntervalDuration_s', 'FinalIntervalCalls', 'MeanFinalVelocity_m_s'});
    newRow.Deceleration_m_s2 = decel_m_s2;
    Tprofiles = [Tprofiles; newRow]; %#ok<AGROW>
end

fprintf('\nExample deceleration profiles for k_r = %g over %.2f--%.2f m:\n', ...
    krExample, d0Scenario_m, dStop_m);
disp(Tprofiles);

%% Plot
fig = figure('Color', 'w', 'Position', [160 160 1000 1000]);
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
hold on; box on; grid on; grid minor;

xPos = 1:numel(krValues);
bw = 0.22;
bar(xPos - bw, 1000 * Tref.TaMin_vr0_s, 0.22, ...
    'FaceColor', [0.30 0.60 0.85], ...
    'EdgeColor', 'none', ...
    'DisplayName', '$T_a$');
bar(xPos, 1000 * Tref.TbMin_vr0_s, 0.22, ...
    'FaceColor', [0.85 0.45 0.25], ...
    'EdgeColor', 'none', ...
    'DisplayName', '$T_b$');
bar(xPos + bw, 1000 * Tref.IpiMin_vr0_s, 0.22, ...
    'FaceColor', [0.45 0.45 0.45], ...
    'EdgeColor', 'none', ...
    'DisplayName', '$\mathrm{IPI}_{\min}$');

set(gca, 'XTick', xPos, 'XTickLabel', compose('%g', krValues));
xlabel('Responsivity coefficient, $k_r$', 'Interpreter', 'latex');
ylabel('Time (ms)', 'Interpreter', 'latex');
title({'\textbf{A. Minimum timing components at $d_{\mathrm{stop}}$}', ...
    '$v_r=0$ reference used to define the theoretical ceiling'}, ...
    'Interpreter', 'latex');
legend('Location', 'best', 'Interpreter', 'latex', 'FontSize', 9, 'Box', 'on');

yyaxis right
plot(xPos, Tref.CrMax_vr0_Hz, 'ko-', ...
    'LineWidth', 1.2, ...
    'MarkerFaceColor', 'w', ...
    'DisplayName', '$C_{r,\max}$');
ylabel('Call rate (Hz)', 'Interpreter', 'latex');
yyaxis left
axis square

nexttile;
hold on; box on; grid on; grid minor;
rowsRef = Tterminal.kr == krValues(1);
plot(Tterminal.vr_m_s(rowsRef), 1000 * Tterminal.FinalIntervalDuration_s(rowsRef), ':', ...
    'Color', [0.15 0.15 0.15], ...
    'LineWidth', 2.2, ...
    'DisplayName', 'shared across $k_r$');
xlabel('Radial closing velocity, $v_r$ (m s$^{-1}$)', 'Interpreter', 'latex');
ylabel(sprintf('Time from %.0f\\%% $C_{r,\\max}$ to $d_{\\mathrm{stop}}$ (ms)', 100 * finalFraction), ...
    'Interpreter', 'latex');
title({'\textbf{B. Velocity controls time spent in the final spatial interval}', ...
    '$\vphantom{v_r=0}$'}, ...
    'Interpreter', 'latex');

text(0.97, 0.94, 'The time in terminal zone is independent of $k_r$', ...
    'Units', 'normalized', ...
    'Interpreter', 'latex', ...
    'FontSize', 10, ...
    'Color', 'k', ...
    'HorizontalAlignment', 'right', ...
    'VerticalAlignment', 'top');
axis square

nexttile;
hold on; box on; grid on; grid minor;
for i = 1:numel(krValues)
    kr = krValues(i);
    rows = Tterminal.kr == kr;
    plot(Tterminal.vr_m_s(rows), Tterminal.NCallsFinalInterval(rows), '-', ...
        'Color', colorsKr(i,:), ...
        'LineWidth', 1.8, ...
        'DisplayName', sprintf('$k_r=%g$', kr));
end
xlabel('Radial closing velocity, $v_r$ (m s$^{-1}$)', 'Interpreter', 'latex');
ylabel(sprintf('Calls from %.0f\\%% $C_{r,\\max}$ to $d_{\\mathrm{stop}}$', 100 * finalFraction), ...
    'Interpreter', 'latex');
title({'\textbf{C. Responsivity controls call density within that interval}', ...
    '$\vphantom{v_r=0}$'}, ...
    'Interpreter', 'latex');
legend('Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 9, 'Box', 'on');
axis square

nexttile;
hold on; box on; grid on; grid minor;
set(gca, 'XDir', 'reverse');
xlim([dStop_m d0Scenario_m]);
ylim([0 3.4]);

xline(dFractionRef_m, ':', 'Color', [0.35 0.35 0.35], 'LineWidth', 1.2, ...
    'HandleVisibility', 'off');
xline(dStop_m, '--', 'Color', [0.45 0.45 0.45], 'LineWidth', 1.2, ...
    'HandleVisibility', 'off');
patch([dStop_m dFractionRef_m dFractionRef_m dStop_m], [0 0 3.4 3.4], ...
    [0.85 0.88 0.93], 'EdgeColor', 'none', 'FaceAlpha', 0.35, ...
    'HandleVisibility', 'off');

for i = 1:height(Tprofiles)
    vStop = Tprofiles.VelocityAtStop_m_s(i);
    frac = (d0Scenario_m - dPath_m) / (d0Scenario_m - dStop_m);
    frac = min(max(frac, 0), 1);
    vPath_m_s = vStart_m_s + (vStop - vStart_m_s) .* frac;

    plot(dPath_m, vPath_m_s, profileStyles{i}, ...
        'Color', profileColors(i,:), ...
        'LineWidth', 2.2, ...
        'DisplayName', char(Tprofiles.Profile(i)));

    % With the reversed distance axis, use a slightly larger x-value so the
    % timing annotations stay inside the visible terminal-interval region.
    xText = min(d0Scenario_m - 0.12, dFractionRef_m + 0.16);
    yBase = interp1(dPath_m, vPath_m_s, dFractionRef_m, 'linear', 'extrap');
    yOffsets = [0.25, 0.25, 0.33, -0.14];
    yText = yBase + yOffsets(i);
    txt = sprintf('%.0f ms; %.1f calls', ...
        1000 * Tprofiles.FinalIntervalDuration_s(i), Tprofiles.FinalIntervalCalls(i));
    text(xText, yText, txt, ...
        'Color', profileColors(i,:), ...
        'FontSize', 8, ...
        'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'middle');
end

legendX = 1.08;
legendYTop = 0.92;
legendLineStep = 0.12;
for i = 1:height(Tprofiles)
    if Tprofiles.Deceleration_m_s2(i) == 0
        legendTxt = 'a = 0';
    else
        legendTxt = sprintf('a = -%.2f m s^{-2}', Tprofiles.Deceleration_m_s2(i));
    end
    text(legendX, legendYTop - (i-1) * legendLineStep, legendTxt, ...
        'Color', profileColors(i,:), ...
        'FontSize', 8, ...
        'Interpreter', 'tex', ...
        'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'middle');
end

xlabel('Anchor distance, $d$ (m)', 'Interpreter', 'latex');
ylabel('Radial closing velocity, $v_r$ (m s$^{-1}$)', 'Interpreter', 'latex');
title({'\textbf{D. Deceleration lengthens the terminal buzz}', ...
    sprintf('$k_r=%g$, final interval %.2f--%.2f m', krExample, dFractionRef_m, dStop_m)}, ...
    'Interpreter', 'latex');
axis square

formatLatex(fig, "full-square");

if savePlots
    exportPaperFigure(fig, fullfile(outDir, 'explore_shortest_Tb_at_dstop'));
end
