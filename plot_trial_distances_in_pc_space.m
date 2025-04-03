function plot_trial_distances_in_pc_space(session_id, window, sessions, scores_by_port, pcs_to_viz)
    sbp25_sm = movmean(scores_by_port.sbp25{session_id}(:, pcs_to_viz), window);
    sbp50_sm = movmean(scores_by_port.sbp50{session_id}(:, pcs_to_viz), window);  
    sbp75_sm = movmean(scores_by_port.sbp75{session_id}(:, pcs_to_viz), window);
    
    sbp25_cds = cumsum(sqrt(sum(diff(sbp25_sm) .^ 2, 2))); 
    sbp50_cds = cumsum(sqrt(sum(diff(sbp50_sm) .^ 2, 2)));
    sbp75_cds = cumsum(sqrt(sum(diff(sbp75_sm) .^ 2, 2)));

    t = tiledlayout(1, 3);

    % divide xline by 3 b/c ESP distributed over 3 ports
    nexttile; 
    plot(sbp25_cds);
    hold on
    plot([1 size(sbp25_cds, 1)], [sbp25_cds(1), sbp25_cds(end)]);
    xline(sessions.esp(session_id)/3, '--'); 
    title("P25");
    hold off

    nexttile; 
    plot(sbp50_cds); 
    hold on
    plot([1 size(sbp50_cds, 1)], [sbp50_cds(1), sbp50_cds(end)]);
    xline(sessions.esp(session_id)/3, '--');
    title("P50");
    hold off

    nexttile; 
    plot(sbp75_cds); 
    hold on
    plot([1 size(sbp75_cds, 1)], [sbp75_cds(1), sbp75_cds(end)]);
    xline(sessions.esp(session_id)/3, '--');
    title("P75");
    hold off

    title(t, "Trial-Trial Distances (session " + sprintf('%.1s', scores_by_port.id{session_id}) + ")"); 
    xlabel(t, "Trial Number");
    ylabel(t, "Cumulative Distance (L2)"); 
