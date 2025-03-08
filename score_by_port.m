function [sbp25, sbp50, sbp75] = score_by_port(i, sessions, top3_pcs)
    p25 = sessions.p25(i); 
    p50 = sessions.p50(i); 
    p75 = sessions.p75(i);

    ports_entered = top3_pcs.port_entries{i};
    pca_scores = top3_pcs.scores{i};
    
    sbp25 = pca_scores(ports_entered == p25, :); 
    sbp50 = pca_scores(ports_entered == p50, :); 
    sbp75 = pca_scores(ports_entered == p75, :); 
end