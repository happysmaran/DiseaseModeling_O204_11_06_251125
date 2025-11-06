function [mu_peak, sigma_peak, all_ts] = run_experiments_agent(params, agents_template, country_agent_ids, G_country, ban_time)
% run_experiments_agent Run n_runs stochastic simulations and return stats.
% Returns mean and std of combined peak infecteds across runs, and optionally
% sample time series (first run) for inspection.

n_runs = params.n_runs;
peaks = zeros(n_runs,1);
all_ts = cell(n_runs,1);
for r = 1:n_runs
    % seed for this run for reproducibility
    if strcmpi(params.seed_mode,'fixed')
        seed = r + 1000;
    else
        seed = 'shuffle';
    end
    [peak, ts] = run_experiment_single(params, agents_template, country_agent_ids, G_country, ban_time, seed);
    peaks(r) = peak;
    all_ts{r} = ts;
end

mu_peak = mean(peaks);
sigma_peak = std(peaks);
% return aggregated ts for first run for convenience
all_ts = all_ts{1};
end