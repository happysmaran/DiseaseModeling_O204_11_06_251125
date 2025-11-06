function [peak_combined, ts] = run_experiment_single(params, agents_template, country_agent_ids, G_country, ban_time, rng_seed)
% run_experiment_single Run a single stochastic agent-based simulation.
% rng_seed can be 'shuffle' or numeric.

if isnumeric(rng_seed)
    rng(rng_seed);
else
    rng('shuffle');
end

% clone agents template
agents = agents_template; % struct array (by value copy)

Nc = params.N_countries;
t_max = params.t_max;
dt = params.dt;
nsteps = floor(t_max/dt) + 1;
days = 0:dt:t_max;

% Travel matrix baseline
M_pre = params.M;

% book-keeping time series
I_hist = zeros(nsteps, Nc);
S_hist = zeros(nsteps, Nc);
R_hist = zeros(nsteps, Nc);

% helper: mapping from agent id to index in agents array (it's simply id)
N_total = numel(agents);

% precompute neighbors per agent from all country graphs (static)
neighbors = cell(N_total,1);
for c = 1:Nc
    G = G_country{c};
    if isempty(G.Edges) && numnodes(G)==0
        continue;
    end
    if numnodes(G) > 0
        nodes = G.Nodes.Name; % these are agent ids as numbers but stored as doubles
        for ni = 1:numel(nodes)
            aid = str2double(nodes{ni});
            nb = neighbors_in_graph(G, aid);
            neighbors{aid} = nb;
        end
    end
end

% initial counts
for t = 1:nsteps
    % record
    countsS = sum([agents.state] == 0);
    countsI = sum([agents.state] == 1);
    countsR = sum([agents.state] == 2);
    % even though counts per country are more useful:
    for c = 1:Nc
        ids = country_agent_ids{c};
        S_hist(t,c) = sum([agents(ids).state] == 0);
        I_hist(t,c) = sum([agents(ids).state] == 1);
        R_hist(t,c) = sum([agents(ids).state] == 2);
    end

    if t == nsteps
        break;
    end

    current_day = days(t);

    % Travel decision for each agent (based on pre-ban or banned matrix)
    if current_day >= ban_time
        M = zeros(Nc);
    else
        M = M_pre;
    end

    % For each agent: decide travel (probability row-sum) and destination
    for a = 1:N_total
        cur_loc = agents(a).loc;
        row = M(cur_loc, :);
        out_prob = sum(row);
        if out_prob > 0 && rand() < out_prob
            % choose destination according to normalized row
            probs = row;
            probs(cur_loc) = 0;
            probs = probs / sum(probs);
            dest = randsample(1:Nc,1,true,probs);
            agents(a).loc = dest;
        else
            % stay
        end
    end

    % Transmission: for each infected agent, attempt to infect susceptible neighbors
    p_trans = params.beta_per_edge;
    for a = 1:N_total
        if agents(a).state ~= 1
            continue;
        end
        % infected agent a can infect neighbors that are currently in same location
        nb = neighbors{a};
        if isempty(nb), continue; end
        % keep neighbors that are currently colocated and susceptible
        colocated_nb = nb([agents(nb).loc] == agents(a).loc);
        if isempty(colocated_nb), continue; end
        for nbid = colocated_nb
            if agents(nbid).state == 0 % susceptible
                if rand() < p_trans
                    agents(nbid).state = 1; % newly infected
                    agents(nbid).days_infected = 0;
                end
            end
        end
    end

    % Recovery and days_infected increment
    for a = 1:N_total
        if agents(a).state == 1
            agents(a).days_infected = agents(a).days_infected + dt;
            if rand() < params.gamma
                agents(a).state = 2; % recovered
            end
        end
    end
end

I_comb = sum(I_hist,2);
peak_combined = max(I_comb);

ts.days = days;
ts.S = S_hist;
ts.I = I_hist;
ts.R = R_hist;
ts.I_comb = I_comb;
end

function nb = neighbors_in_graph(G, aid)
% G.Nodes.Name are strings (when constructed using numeric ids)
nodelist = G.Nodes.Name;
idx = find(strcmp(nodelist, num2str(aid)));
if isempty(idx)
    nb = [];
    return;
end
nei = neighbors(G, idx);
% map neighbor node indices to agent ids (stored in G.Nodes.Name)
nb = zeros(numel(nei),1);
for k = 1:numel(nei)
    nb(k) = str2double(G.Nodes.Name{nei(k)});
end
end
