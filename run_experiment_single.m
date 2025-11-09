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
        nodes = G.Nodes.Name; % cellstr or string of agent ids
        for ni = 1:numel(nodes)
            aid = str2double(nodes{ni});
            nb = neighbors_in_graph(G, aid);
            neighbors{aid} = nb;
        end
    end
end

% % --- DIAGNOSTIC BLOCK (paste after neighbors built) ---
% fprintf('DIAG: checking neighbors and agent locs...\n');
% 
% % basic agent / country checks
% Nc = params.N_countries;
% N_total = numel(agents);
% locs = [agents.loc];
% if any(locs < 1) || any(locs > Nc) || any(mod(locs,1)~=0)
%     warning('DIAG: some agents have invalid loc values (outside 1..Nc or non-integer)');
% end
% 
% % neighbor validity checks
% invalid_count = 0;
% nan_count = 0;
% out_of_range = 0;
% deg = zeros(N_total,1);
% for a = 1:N_total
%     nb = neighbors{a};
%     deg(a) = numel(nb);
%     if isempty(nb), continue; end
%     if any(isnan(nb))
%         nan_count = nan_count + 1;
%     end
%     if any(nb < 1) || any(nb > N_total) || any(mod(nb,1)~=0)
%         out_of_range = out_of_range + 1;
%     end
%     % also check that neighbor ids are unique and integers
%     if any(diff(sort(nb))==0)
%         % duplicate neighbor ids (not fatal but suspicious)
%     end
% end
% 
% fprintf('DIAG: N_total=%d, Nc=%d\n', N_total, Nc);
% fprintf('DIAG: degree mean=%.3f median=%.3f max=%d zero_deg=%d\n', mean(deg), median(deg), max(deg), sum(deg==0));
% fprintf('DIAG: neighbor lists with NaNs = %d, out_of_range lists = %d\n', nan_count, out_of_range);
% 
% % show examples for first few infected (or first few agents)
% example_ids = unique([params.I0_idx, find([agents.state]==1)]);
% example_ids = example_ids(1:min(10,length(example_ids)));
% for id = example_ids
%     if id > N_total, continue; end
%     nb = neighbors{id};
%     if isempty(nb)
%         fprintf('DIAG: agent %d (state=%d, loc=%d) has ZERO neighbors\n', id, agents(id).state, agents(id).loc);
%     else
%         % show first few neighbor ids and their locs and states
%         nshow = min(8,numel(nb));
%         nb_small = nb(1:nshow);
%         nb_locs = [agents(nb_small).loc];
%         nb_states = [agents(nb_small).state];
%         fprintf('DIAG: agent %d (loc=%d) neighbors [%s] locs [%s] states [%s]\n', ...
%             id, agents(id).loc, num2str(nb_small), num2str(nb_locs), num2str(nb_states));
%     end
% end
% 
% % Check if neighbor ids are local (i.e., only span 1..max_local)
% all_node_ids_in_graphs = [];
% for c=1:Nc
%     G = G_country{c};
%     if numnodes(G)>0
%         node_names = G.Nodes.Name;
%         % try converting to numbers
%         nums = str2double(node_names);
%         all_node_ids_in_graphs = [all_node_ids_in_graphs; nums(:)];
%     end
% end
% fprintf('DIAG: node id range seen in graphs: min=%g max=%g (NaNs=%d)\n', min(all_node_ids_in_graphs), max(all_node_ids_in_graphs), sum(isnan(all_node_ids_in_graphs)));
% % --- END DIAGNOSTIC BLOCK ---



% main time loop
for t = 1:nsteps
    % record per-country counts (based on current agents)
    % record per-country counts (based on current agents' loc)
    locs = [agents.loc];                 % 1 x N_total vector of current locations
    for c = 1:Nc
        ids = find(locs == c);           % current agent ids in country c
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
    % NOTE: travel changes loc in-place; travel is independent from infection decision in this model.
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

    % SYNCHRONOUS INFECTION & RECOVERY
    old_state = [agents.state];  % snapshot at start of step
    will_be_infected = false(N_total,1); % mark new infections (to apply after checking all infectious)
    will_recover = false(N_total,1);     % mark recoveries (to apply after decisions)
    p_trans = params.beta_per_edge;

    % Infection decisions: infected (in old_state==1) attempt to infect suscept in same loc
    for a = 1:N_total
        if old_state(a) ~= 1
            continue;
        end
        nb = neighbors{a};
        if isempty(nb), continue; end
        % only consider neighbors that are colocated and susceptible *at start of step*
        colocated_nb = nb([agents(nb).loc] == agents(a).loc);
        if isempty(colocated_nb), continue; end
        for nbid = colocated_nb'
            if old_state(nbid) == 0 % susceptible at start of step
                if rand() < p_trans
                    will_be_infected(nbid) = true;
                end
            end
        end
    end

    % Recovery decisions: infected at start of step may recover during this step
    for a = 1:N_total
        if old_state(a) == 1
            % sample recovery using params.gamma (probability per day)
            if rand() < params.gamma
                will_recover(a) = true;
            end
        end
    end

    % Apply updates synchronously:
    %  - those marked to recover become recovered
    %  - those marked to be infected become infected unless they were already recovered (defensive)
    for a = 1:N_total
        if will_recover(a)
            agents(a).state = 2; % recovered
            agents(a).days_infected = 0;
        end
    end

    % Apply infections AFTER recoveries so a person who recovered this step doesn't get re-infected
    for a = 1:N_total
        if will_be_infected(a)
            if agents(a).state == 0  % only infect susceptibles
                agents(a).state = 1;
                agents(a).days_infected = 0;
            end
        end
    end

    % increment days_infected for those still infected
    for a = 1:N_total
        if agents(a).state == 1
            agents(a).days_infected = agents(a).days_infected + dt;
        end
    end
    % END synchronous update
end

infected_countries = find(max(I_hist) > 0);
%fprintf('Infected countries: %d of %d\n', numel(infected_countries), Nc);
%disp(infected_countries);

I_comb = sum(I_hist,2);
peak_combined = max(I_comb);

ts.days = days;
ts.S = S_hist;
ts.I = I_hist;
ts.R = R_hist;
ts.I_comb = I_comb;

%total_infected = sum(ts.R(end,:) + ts.I(end,:));
%fprintf('Total ever infected = %d of %d\n', total_infected, N_total);

end