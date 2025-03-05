function [sbp25, sbp50, sbp75] = score_by_port(i, sessions, top3_pcs)
    esp = sessions.esp(i);
    p25 = sessions.p25(i); 
    p50 = sessions.p50(i); 
    p75 = sessions.p75(i); 
    ports_entered = top3_pcs.port_entries{i};

    trials25 = find(ports_entered(1 : esp-1) == p25); 
    trials25 = [trials25; (find(ports_entered(esp : end) == p75))];

    trials50 = find(ports_entered == p50);

    trials75 = find(ports_entered(1 : esp-1) == p75);
    trials75 = [trials75; (find(ports_entered(esp : end) == p25))];
    
    pca_scores = top3_pcs.scores{i}; 
    sbp25 = pca_scores(trials25, :); 
    sbp50 = pca_scores(trials50, :); 
    sbp75 = pca_scores(trials75, :); 
    clearvars -except sbp25 sbp50 sbp75
end