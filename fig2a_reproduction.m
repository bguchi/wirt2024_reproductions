%% Preliminaries
% Data Selection

% extract session .mat names
listing = dir("3portdata"); 
listing = struct2table(listing);
session_ids = listing.name(~listing.isdir); 

% select only "good" sessions (per Dr. Hyman)
good_sessions = readmatrix("good_sessions.csv"); 
session_ids = string(session_ids(good_sessions));

% extract recorded session ESPs
load("SwitchTrials.mat", "SwTrial"); 
session_esps = SwTrial(good_sessions, 2); 
clear listing good_sessions SwTrial

% Port Inference
% Per |hyman2017_reproductions/group.mlx| port inference routine: 

sz = [size(session_ids, 1) 5];
vartypes = ["string" repmat("double", [1 4])]; 
varnames = ["id", "esp", "p25", "p50", "p75"]; 
sessions = table('Size', sz, 'VariableTypes', vartypes, 'VariableNames', varnames); 

for i = 1 : height(sessions)
    load("3portdata/" + session_ids(i), "Event_timestamps"); 
    ports = port_ids(Event_timestamps);

    p25 = find(ismember(ports, "25%"));
    p50 = find(ismember(ports, "50%"));
    p75 = find(ismember(ports, "75%"));
    sessions(i, :) = {session_ids(i), session_esps(i), p25, p50, p75};
    clear i Event_timestamps ports p25 p50 p75
end 
clear sz vartypes varnames session_ids session_esps

%% PCA Reproduction: Fig. 2A
% Running PCA on iFR data filtered for NP events and averaged over event durations: 

pcs = cell(height(sessions), 5); 
for i = 1 : height(sessions) 
    [filtered_and_avgd_iFRs, np_ports] = filter_and_avg_iFRs(sessions.id(i)); 
    [coeff, score, ~, ~, explained, ~] = pca(filtered_and_avgd_iFRs); 
    pcs{i, 1} = sessions.id(i); 
    pcs{i, 2} = np_ports; 
    pcs{i, 3} = coeff; 
    pcs{i, 4} = score; 
    pcs{i, 5} = explained; 
    clear i filtered_and_avgd_iFRs np_ports coeff score explained
end
pcs = cell2table(pcs, 'VariableNames', ["id", "port_entries", "coeffs", "scores", "var_exp"])

% Comparing average variances explained by each top 3 PC to reported values: 
explained = zeros(height(sessions), 3); 
for i = 1 : height(sessions)
    var_exp = pcs.var_exp{i}'; 
    explained(i, :) = var_exp(1:3); 
end

disp("Average var. explained by PC1: " + mean(explained(:, 1)) + "; reported 9.72%")
disp("Average var. explained by PC2: " + mean(explained(:, 2)) + "; reported 5.93%")
disp("Average var. explained by PC3: " + mean(explained(:, 3)) + "; reported 5.17%")

total = sum(mean(explained, 1)); 
disp("Total var. explained by top 3 PCs: " + total + "; reported 20.8%")
clear i explained total 

% Faceting PCA scores by port contingency for visualization:
scores_by_port = cell(height(sessions), 4); 
for i = 1 : height(sessions)
    scores_by_port{i, 1} = sessions.id(i); 
    [sbp25, sbp50, sbp75] = score_by_port(i, sessions, pcs);
    scores_by_port{i, 2} = sbp25; 
    scores_by_port{i, 3} = sbp50; 
    scores_by_port{i, 4} = sbp75;
    clear i sbp25 sbp50 sbp75 
end
scores_by_port = cell2table(scores_by_port, 'VariableNames', ["id", "sbp25", "sbp50", "sbp75"])

% Plotting!
session_id = 7;
window = 15;
pcs_to_viz = 1:3;
ports_to_viz = ["p25" "p50" "p75"]; 

plot_smooth_trajectories(session_id, window, scores_by_port, pcs_to_viz, ports_to_viz)
