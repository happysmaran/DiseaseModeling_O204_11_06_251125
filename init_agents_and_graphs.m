function [agents, country_agent_ids, G_country] = init_agents_and_graphs(params)
% init_agents_and_graphs
% -------------------------------------------------------------
% Create agent list and per-country static contact graphs.
% Each agent has:
%   id              - unique integer ID
%   home            - home country (1..Nc)
%   loc             - current country
%   state           - 0=S, 1=I, 2=R
%   days_infected   - time since infection
%
% Also builds per-country contact graphs (Erdős–Rényi).
% -------------------------------------------------------------

Nc = params.N_countries;
pop = params.population;
N_total = sum(pop);

ids = num2cell((1:N_total)');  % column cell vector (avoids struct dimension error)
home = num2cell(zeros(N_total,1));
loc = num2cell(zeros(N_total,1));
state = num2cell(zeros(N_total,1));
days_infected = num2cell(zeros(N_total,1));

agents = struct('id', ids, ...
                'home', home, ...
                'loc', loc, ...
                'state', state, ...
                'days_infected', days_infected);

% Assign agents to countries
country_agent_ids = cell(Nc,1);
idx = 1;
for c = 1:Nc
    n = pop(c);
    ids_c = idx:(idx+n-1);
    country_agent_ids{c} = ids_c;
    idx = idx + n;
    for j = ids_c
        agents(j).home = c;
        agents(j).loc = c;  % everyone starts at home
        agents(j).state = 0; % susceptible initially
        agents(j).days_infected = 0;
    end
end

% Seed initial infected agent(s)
I0 = params.I0_idx;
if isscalar(I0)
    if I0 < 1 || I0 > N_total
        I0 = 1;
    end
    agents(I0).state = 1;
    agents(I0).days_infected = 0;
else
    for ii = I0
        if ii >= 1 && ii <= N_total
            agents(ii).state = 1;
            agents(ii).days_infected = 0;
        end
    end
end

% Build per-country contact graphs (Erdős–Rényi)
G_country = cell(Nc,1);
for c = 1:Nc
    ids_c = country_agent_ids{c};
    nc = numel(ids_c);
    if nc <= 1
        G_country{c} = graph();
        continue;

    end

    p = params.p_edge;  % edge probability
    A = rand(nc) < p;

    A = triu(A,1);      % upper triangle only
    A = A + A';         % make symmetric (undirected)
    names = string(ids_c);
    G_country{c} = graph(A, names);  % node names are agent IDs
end
end