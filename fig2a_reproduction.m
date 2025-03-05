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
% Per |hyman2017_reproductions/group.mlx| routine: 

sz = [size(session_ids, 1) 5];
vartypes = ["string" repmat("double", [1 4])]; 
varnames = ["id", "esp", "p25", "p50", "p75"]; 
sessions = table('Size', sz, 'VariableTypes', vartypes, 'VariableNames', varnames); 

for i = 1 : height(sessions)
    load("3portdata/" + session_ids(i), "Event_timestamps"); 
    ports = port_ids(Event_timestamps, 10);

    p25 = find(ismember(ports, "25%"));
    p50 = find(ismember(ports, "50%"));
    p75 = find(ismember(ports, "75%"));
    sessions(i, :) = {session_ids(i), session_esps(i), p25, p50, p75};
    clear i Event_timestamps ports p25 p50 p75
end 
clear sz vartypes varnames session_ids session_esps
%% PCA Reproduction: Fig. 2A
% Running PCA on iFR data filtered for NP events and averaged over event durations: 

top3_pcs = cell(height(sessions), 4); 
for i = 1 : height(sessions) 
    [filtered_and_avgd_iFRs, np_ports] = filter_and_avg_iFRs(sessions.id(i)); 
    [coeff, score, ~, ~, explained, ~] = pca(filtered_and_avgd_iFRs); 
    top3_pcs{i, 1} = sessions.id(i); 
    top3_pcs{i, 2} = np_ports; 
    top3_pcs{i, 3} = coeff(:, 1:3); 
    top3_pcs{i, 4} = score(:, 1:3); 
    top3_pcs{i, 5} = explained(1:3); 
    clear i filtered_and_avgd_iFRs np_ports coeff score explained
end
top3_pcs = cell2table(top3_pcs, 'VariableNames', ["id", "port_entries", "coeffs", "scores", "var_exp"])
%% 
% Comparing average variances explained by each PC to reported values: 

explained = zeros(height(sessions), 3); 
for i = 1 : height(sessions)
    explained(i, :) = top3_pcs.var_exp{i}'; 
end

disp("Average var. explained by PC1: " + mean(explained(:, 1)) + "; reported 9.72%")
disp("Average var. explained by PC2: " + mean(explained(:, 2)) + "; reported 5.93%")
disp("Average var. explained by PC3: " + mean(explained(:, 3)) + "; reported 5.17%")

total = sum(mean(explained, 1)); 
disp("Total var. explained by top 3 PCs: " + total + "; reported 20.8%")
clear i explained total 
%% 
% Faceting PCA scores by port contingency for visualization:

scores_by_port = cell(height(sessions), 4); 
for i = 1 : height(sessions)
    scores_by_port{i, 1} = sessions.id(i); 
    [sbp25, sbp50, sbp75] = score_by_port(i, sessions, top3_pcs);
    scores_by_port{i, 2} = sbp25; 
    scores_by_port{i, 3} = sbp50; 
    scores_by_port{i, 4} = sbp75; 
    clear i sbp25 sbp50 sbp75
end
scores_by_port = cell2table(scores_by_port, 'VariableNames', ["id", "sbp25", "sbp50", "sbp75"])
%% 
% Visualizing trajectories of neural activity in normalized top 3 PC space, 
% faceted by port contingency: 

% visualize 3d plot for session 2 ('b', index 1)
session_id = 1; 
sbp25 = normalize(scores_by_port.sbp25{session_id});
sbp50 = normalize(scores_by_port.sbp50{session_id}); 
sbp75 = normalize(scores_by_port.sbp75{session_id});

plot3(sbp25(:, 1), sbp25(:, 2), sbp25(:, 3), ...
      sbp50(:, 1), sbp50(:, 2), sbp50(:, 3), ...
      sbp75(:, 1), sbp75(:, 2), sbp75(:, 3))
xlabel('PC1')
ylabel('PC2')
zlabel('PC3')
title("session " + sprintf('%.1s', sessions.id(session_id)))
legend(["25/75" "50/50", "75/25"])
clear session_id sbp25 sbp50 sbp75