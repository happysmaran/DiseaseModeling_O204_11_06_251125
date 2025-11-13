function nb = neighbors_in_graph(G, aid)
% Return neighbor agent ids (numeric vector) for agent 'aid' in graph G.
% Robust to G.Nodes.Name being cellstr or string array.

nodelist = G.Nodes.Name;
% Convert to cellstr for uniform strcmp usage
if isstring(nodelist)
    nodelist_cs = cellstr(nodelist);
elseif iscell(nodelist)
    nodelist_cs = nodelist;
else
    nodelist_cs = cellstr(string(nodelist));
end

idx = find(strcmp(nodelist_cs, num2str(aid)));
if isempty(idx)
    nb = [];
    return;
end
nei_idx = neighbors(G, idx); % returns indices of neighbor nodes
nb = zeros(numel(nei_idx),1);
for k = 1:numel(nei_idx)
    nb(k) = str2double(nodelist_cs{nei_idx(k)});
end
end