%% find the Tb* and time to contact

close all
% Define parameters
type = "fm";
motile = 0;
makeAudio = 0; % do not synthesize call sequence
kr = [3:7];                % Time delay factor
initial_velocity = [3:7];            % Prey velocity towards bat (m/s)
if type == "fm"
    target_distance = 10;            % Initial target distance (meters)
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
%% Create summary table
sum_table = struct();
for i = 1: height(T)
    sum_table(i).kr = T.Kr{i};
    sum_table(i).vel= T.initial_velocity{i};
    sum_table(i).Tb_idx = T.Tb_idx{i};
    ttc = T.time_to_contact{i}*1000; % in ms
    sum_table(i).time_to_contact = ttc(sum_table(i).Tb_idx);
    dtt =  T.new_target_distance{i};
    sum_table(i).distance_rmax = dtt(sum_table(i).Tb_idx);
    sum_table(i).Tb_prime = abs(T.Tb_prime{i}*1000); % in ms
    sum_table(i).Tb = T.Tb{i}*1000; % in ms
end
%%
% Extract values from struct array
kr = [sum_table.kr]';
vel = [sum_table.vel]';
Tb_prime = [sum_table.Tb_prime]';
Ttc = [sum_table.time_to_contact]';
Tb = [sum_table.Tb]';
Dtt = [sum_table.distance_rmax]';
% Define fine grid
krq = linspace(min(kr), max(kr), 500);
velq = linspace(min(vel), max(vel), 500);
[KQ, VQ] = meshgrid(krq, velq);

% Interpolate each quantity using cubic interpolation
Tb_prime_interp = griddata(kr, vel, Tb_prime, KQ, VQ, 'cubic');
Ttc_interp = griddata(kr, vel, Ttc, KQ, VQ, 'cubic');
Tb_interp = griddata(kr, vel, Tb, KQ, VQ, 'cubic');
Dtt_interp = griddata(kr, vel, Dtt, KQ, VQ, 'cubic');
% Plot all three as contour plots
figure('Position', [100 100 1400 400]);

subplot(2,2,1)
contourf(KQ, VQ, Tb_prime_interp, 20, 'EdgeColor','none');
cb = colorbar;
cb.TickLabelInterpreter = 'latex';
cb.Label.String = 'ms';
cb.Label.Interpreter = 'latex';
xlabel('$k_r$');
ylabel('Velocity (m/s)');
title('$T_{b^*}$ (ms)');
axis square
formatLatex(gca);

subplot(2,2,2)
contourf(KQ, VQ, Ttc_interp, 20, 'EdgeColor','none');
cb = colorbar;
cb.TickLabelInterpreter = 'latex';
cb.Label.String = 'ms';
cb.Label.Interpreter = 'latex';
xlabel('$k_r$');
ylabel('Velocity (m/s)');
title("Time to Contact (ms)");
axis square
formatLatex(gca);

subplot(2,2,3)
contourf(KQ, VQ, Tb_interp, 20, 'EdgeColor','none');
cb = colorbar;
cb.TickLabelInterpreter = 'latex';
cb.Label.String = 'ms';
cb.Label.Interpreter = 'latex';
xlabel('$k_r$');
ylabel('Velocity (m/s)');
title("$T_b$ at $\mathcal{R}_{max}$ (ms)");
axis square
formatLatex(gca);

subplot(2,2,4)
contourf(KQ, VQ, Dtt_interp, 20, 'EdgeColor','none');
cb = colorbar;
cb.TickLabelInterpreter = 'latex';
cb.Label.String = 'm';
cb.Label.Interpreter = 'latex';
xlabel('$k_r$');
ylabel('Velocity (m/s)');
title("Dist. to Target at $\mathcal{R}_{max}$ (m)");
axis square
formatLatex(gca);

sgtitle("Spatio-Temporal Biosonar Parameters", 'Interpreter', 'latex', ...
                              'FontSize', 18, ...
                              'FontWeight', 'bold');


%%
function T = resultStructToRow(result)
fn = fieldnames(result);
values = cell(1, numel(fn));

for i = 1:numel(fn)
    values{i} = {result.(fn{i})};  % Wrap each field in a cell
end

T = cell2table(values, 'VariableNames', fn);
end