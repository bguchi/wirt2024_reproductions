function plot_smooth_trajectories(session_id, window, scores_by_port)
    sbp25_sm = movmean(scores_by_port.sbp25{session_id}, window);
    sbp50_sm = movmean(scores_by_port.sbp50{session_id}, window);  
    sbp75_sm = movmean(scores_by_port.sbp75{session_id}, window);
    
    figure
    ax = axes(); 
    limx = max(abs([sbp25_sm(:, 1); sbp50_sm(:, 1); sbp75_sm(:, 1)]));
    limy = max(abs([sbp25_sm(:, 2); sbp50_sm(:, 2); sbp75_sm(:, 2)]));
    limz = max(abs([sbp25_sm(:, 3); sbp50_sm(:, 3); sbp75_sm(:, 3)]));
    
    xlim(ax, [-ceil(limx) ceil(limx)]);
    ylim(ax, [-ceil(limy) ceil(limy)]);
    zlim(ax, [-ceil(limz) ceil(limz)]);
    view(ax, 3)

    sbps = {sbp25_sm sbp50_sm sbp75_sm}; 
    for i = 1 : size(sbps, 2)
        x = sbps{i}(:, 1); 
        y = sbps{i}(:, 2);
        z = sbps{i}(:, 3); 
        z(end) = NaN; 
        c = (1 : size(sbps{i}, 1))'; 
        patch(x, y, z, c, 'EdgeColor','interp','LineWidth', 2, 'LineJoin','round');
    end

    xlabel('PC1'); 
    ylabel('PC2'); 
    zlabel('PC3'); 
    title("session " + sprintf('%.1s', scores_by_port.id{session_id})); 
    cb = colorbar; 
    cb.Label.String = 'Trial Number';
end