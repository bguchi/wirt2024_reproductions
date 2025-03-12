function [sbp25, sbp50, sbp75] = score_by_port(i, sessions, pcs)
    % SCORE_BY_PORT Arranges PC scores by corresponding port entry. 
    %    Args: 
    %        - i: (int) to index SESSIONs
    %        - sessions: (table) as constructed in `Port Inference` section
    %        - pcs: (table) as constructed in `PCA Reproductions` section
    %    Out: 
    %        - [sbp25, ..., sbp75]: PC scores corresponding to e.g., P25
    p25 = sessions.p25(i); 
    p50 = sessions.p50(i); 
    p75 = sessions.p75(i);

    ports_entered = pcs.port_entries{i};
    pca_scores = pcs.scores{i};
    
    sbp25 = pca_scores(ports_entered == p25, :); 
    sbp50 = pca_scores(ports_entered == p50, :); 
    sbp75 = pca_scores(ports_entered == p75, :); 
end