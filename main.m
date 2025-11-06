clearvars; close all; clc;
rng(1);

% Parameters
params.N_countries = 10;
params.N_total = 100;          % total agents
params.beta_per_edge = 0.15;   % transmission probability per contact per day
params.gamma = 0.10;    
       % recovery prob per day
params.p_edge = 0.05;          % prob edge exists in Erdos-Renyi contact graph within country
params.avg_out_travel = 0.02;  % fraction of persons leaving their country per day (row-sum of M)
params.t_max = 160;            % days
params.dt = 1;                 % day step (1)
params.I0_idx = 1;             % initial infected agent id index (or array)
params.n_runs = 20;            % stochastic repeats per ban time for averaging
params.seed_mode = 'auto';     % 'auto' or 'fixed' (if fixed will use run index to seed)
params.verbose = true;

% Distribute populations: equal by default
Nc = params.N_countries;
params.population = floor(params.N_total / Nc) * ones(Nc,1);
% if remainder, add 1 to first few countries
remn = params.N_total - sum(params.population);
params.population(1:remn) = params.population(1:remn) + 1;

% Build baseline travel matrix M (rows sum to avg_out_travel)
base_rate = params.avg_out_travel;
M = (base_rate/(Nc-1)) * (ones(Nc)-eye(Nc));
params.M = M;

% Ban sweep
ban_search.t_min = 0;
ban_search.t_max = 60;
ban_search.t_step = 1; % daily checks
ban_times = ban_search.t_min:ban_search.t_step:ban_search.t_max;

target_reduction = 50; % percent

% Generate static contact graphs & initial agent assignment
[agents_template, country_agent_ids, G_country] = init_agents_and_graphs(params);
% agents_template is a struct array that can be cloned for each run

% Baseline (no ban)
disp('Running baseline (no travel ban)...');
[peak_baseline_mean,peak_baseline_std, baseline_ts] = run_experiments_agent(params, agents_template, country_agent_ids, G_country, Inf);

% Plot baselne ts
% Plot baseline time series
figure('Name','Baseline Time Series');
plot(baseline_ts.days, baseline_ts.I_comb, 'LineWidth', 2);
xlabel('Days'); ylabel('Combined Infected');
title('Baseline Infection Over Time');
grid on;

fprintf('Baseline mean combined peak infected (over %d runs) = %.3f (std %.3f)\n', params.n_runs, peak_baseline_mean, peak_baseline_std);

% Sweep ban times
n_b = numel(ban_times);
mean_peaks = zeros(n_b,1);
std_peaks  = zeros(n_b,1);

if params.verbose, disp('Sweeping ban times...'); end
for k = 1:n_b
    bt = ban_times(k);
    [mu, sigma, ~] = run_experiments_agent(params, agents_template, country_agent_ids, G_country, bt);
    mean_peaks(k) = mu;
    std_peaks(k) = sigma;
    if params.verbose
        fprintf('Ban at %2d d: mean peak = %.3f, std = %.3f\n', bt, mu, sigma);
    end
end

% compute percent reductions using mean peaks
reductions = 100 * (peak_baseline_mean - mean_peaks) / peak_baseline_mean;

% find earliest ban time achieving target reduction
idx = find(reductions >= target_reduction, 1, 'first');

if isempty(idx)
    fprintf('No ban time in [%d,%d] days reached >= %.1f%% mean peak reduction.\n', ban_search.t_min, ban_search.t_max, target_reduction);
else
    fprintf('Earliest ban day achieving >= %.1f%% mean peak reduction: %d days\n', target_reduction, ban_times(idx));
    fprintf('Mean peak then = %.3f (baseline = %.3f)\n', mean_peaks(idx), peak_baseline_mean);
end

% Plot results: reduction vs ban time and baseline vs best
figure('Name','Peak reduction vs Ban time');
errorbar(ban_times, reductions, 0.5*std_peaks./peak_baseline_mean*100, '-o','LineWidth',1.4);
xlabel('Ban time (days)'); ylabel('Mean peak reduction (%)');
yline(target_reduction, ':r', sprintf('Target = %.0f%%', target_reduction));
title('Mean peak reduction vs travel-ban day');
grid on;

    
if ~isempty(idx)
    best_bt = ban_times(idx);
    % simulate example trajectories for baseline and best ban (single representative run)
    [~, ~, ts_base] = run_experiment_single(params, agents_template, country_agent_ids, G_country, Inf, 12345);
    [~, ~, ts_best] = run_experiment_single(params, agents_template, country_agent_ids, G_country, best_bt, 54321);

    figure('Name','Combined infected (example runs)');
    plot(ts_base.days, ts_base.I_comb, 'LineWidth',2); hold on;
    plot(ts_best.days, ts_best.I_comb, '--','LineWidth',2);
    xlabel('Days'); ylabel('Combined infected');
    legend('Baseline (no ban)','Ban at best day','Location','northeast');
    title(sprintf('Example runs: ban at %d days', best_bt));
    grid on;
end

disp('Done.');