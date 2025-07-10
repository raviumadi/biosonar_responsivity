%% generate biological reaction time plots
close all
% Define parameters
type = "fm";
motile = 0;  
makeAudio = 0;
kr = [4,5,6,7];                % Time delay factor
initial_velocity = [3,4,5];            % Prey velocity towards bat (m/s)
if type == "fm"
    target_distance = 20;            % Initial target distance (meters)
    initial_call_duration = 0.005; % Initial call duration (seconds)
    bandwidth = 65e3:-100:35e3; % Bandwidth sweep (Hz) 
elseif type == "cf"
    target_distance = 70;            % Initial target distance (meters)
    initial_call_duration = 0.05; % Initial call duration (seconds)
    bandwidth = 70e3;
end
T = table();
first = true;

for i = 1:length(kr)
    for j = 1:length(initial_velocity)
        try
            if type == "fm"
                result = simulateEcholocation(bandwidth, kr(i), target_distance, initial_velocity(j), initial_call_duration, motile, makeAudio);
            elseif type == "cf"
                result = simulateEcholocationCF(bandwidth, kr(i), target_distance, initial_velocity(j), initial_call_duration);
            end
            t = resultStructToRow(result);
            t.Kr = {kr(i)};
            t.InitialVelocity = {initial_velocity(j)};
            t.TargetDistance = {target_distance};
            
            if first
                T = t;
                first = false;
            else
                T = [T; t];
            end
        catch ME
            warning('Failed at kr=%d, velocity=%d: %s', kr(i), initial_velocity(j), ME.message);
        end
    end
end

%% for CF calls, plot call duration contractio at buzz point
for i = 1 : height(T)
    call_durations = T.new_call_duration{i};
    buzz_point = T.MaxCallRatePoint{i};
    tb_point = T.Tb_idx{i};
    contracted_call_duration(i) = call_durations(buzz_point);
    contraction_at_tb(i) = call_durations(tb_point);
end
pecentrage_call_contraction = (initial_call_duration - contracted_call_duration)./initial_call_duration*100;
pecentrage_call_contraction_tb = (initial_call_duration - contraction_at_tb)./initial_call_duration*100
%% plot tb = kr * ta
% Define color map for kr = 4:7

kr_colors = {
    [51, 102, 204]./255;    % blue
    [255, 102, 0]./255; % orange
    [204, 0, 153]./255; % magenta
    [102, 102, 51]./255; % olive
};

% Extract unique velocities and kr values
velocities = unique(cell2mat(T.InitialVelocity));
kr_values = unique(cell2mat(T.Kr));

figure; clf
tl = tiledlayout(1, length(velocities), 'TileSpacing','compact');

for v = 1:length(velocities)
   nt =  nexttile;
    v_mask = cellfun(@(v_) v_ == velocities(v), T.InitialVelocity);

    for k = 1:length(kr_values)
        k_mask = cellfun(@(k_) k_ == kr_values(k), T.Kr);
        idx = find(v_mask & k_mask);

        if isempty(idx)
            continue
        end

        % Tb = T.Ta{idx} + T.Ta{idx} * T.Kr{idx};
        Tb = T.Ta{idx} * T.Kr{idx};
        c = kr_colors{k};

        plot(nt, Tb, 'LineWidth', 1.5, 'Color', c, 'DisplayName', ['kr = ' num2str(kr_values(k))]);
        hold on
        plot(nt, T.Tb_idx{idx}, T.Tb{idx}, '.', 'MarkerSize', 15, 'Color', c)
        plot(nt, T.OverlapPoint{idx}, Tb(T.OverlapPoint{idx}), 'r.', 'MarkerSize', 15)
        plot(nt, T.MaxCallRatePoint{idx}, Tb(T.MaxCallRatePoint{idx}), 'k.', 'MarkerSize', 15)
    end
    drawnow
    nt.FontWeight = 'bold';
    nt.TickLabelInterpreter = 'latex';
    nt.FontSize = 18;
    set(nt, 'TickLabelInterpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold');
    title(['Initial velocity = ' num2str(velocities(v)) ' m/s'], 'Interpreter','latex', 'FontSize', 16)
    xlabel('Call Number', 'Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold')
    ylabel('$T_b$ (s)', 'Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold')
    grid on
    grid minor
    ylim([-0.05 0.3])
    axis square

end

% title(tl, "Reaction Time Window $T_b = k_r \cdot T_a$", 'Interpreter','latex', 'FontSize', 18)
if type == "fm"
subtitle(tl, 'FM Bats: Grouped by Initial Velocity, Colored by $k_r$,  Buzz Call Rate - 200Hz', 'Interpreter','latex', 'FontSize', 16)
elseif type == "cf"
subtitle(tl, 'CF Bats : $t_{call}$ - 50ms, Target Distance - 60m, Buzz Call Rate - 100Hz', 'Interpreter','latex', 'FontSize', 16)
end


%% plot call rate and velocity relationhip
for i = 1: height(T)
    dist = T.delta_s{i};
    time_diff = T.delta_t{i};
    time_diff = time_diff(2:end);
    vel = abs(dist./time_diff);
    cr = T.call_rate_profile{i};
    % cr = cr(2:end);
    % filter
    % idx = cr<200 ;

    plot(cr, vel, 'k.');
    hold on
    % plot(vel, 'm-.' )

    % pause
end
%%


function T = resultStructToRow(result)
fn = fieldnames(result);
values = cell(1, numel(fn));

for i = 1:numel(fn)
    values{i} = {result.(fn{i})};  % Wrap each field in a cell
end

T = cell2table(values, 'VariableNames', fn);
end