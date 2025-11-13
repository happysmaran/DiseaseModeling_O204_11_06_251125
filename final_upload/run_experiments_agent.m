function [mu_peak, sigma_peak, all_ts, avg_ts] = run_experiments_agent(params, agents_template, country_agent_ids, G_country, ban_time)
% run_experiments_agent Run n_runs stochastic simulations and return stats.
% Returns mean and std of combined peak infecteds across runs=

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


avg_days = zeros(length(all_ts{1}.I), 1);


for i = 1:n_runs
    avg_days = avg_days + sum(all_ts{i}.I, 2);
    
end

avg_days = avg_days ./ n_runs;


avg_ts = avg_days;
end