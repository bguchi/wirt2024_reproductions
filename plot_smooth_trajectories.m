function plot_smooth_trajectories(session_id, window, scores_by_port, pcs_to_viz, ports_to_viz)
    % PLOT_SMOOTH_TRAJECTORIES of PC activity in 3D. 
    %    Args: 
    %        - session_id: (int) to index SCORES_BY_PORT
    %        - window: (int) "smoothing factor" to `movmean`
    %        - scores_by_port: (table) as constructed in `PCA Reproductions`
    %        - pcs_to_viz: (1 x 3 int) PCs to visualize
    %        - ports_to_viz: (1 x 1-3 string) port trajectories to visualize; 
    %          ports must be of form "p{25, 50, 75}"
    %    Out: 
    %        - plot3 figure of smoothed trajectories, color-coded for trial
    %          number
    sbp25_sm = movmean(scores_by_port.sbp25{session_id}(:, pcs_to_viz), window);
    sbp50_sm = movmean(scores_by_port.sbp50{session_id}(:, pcs_to_viz), window);  
    sbp75_sm = movmean(scores_by_port.sbp75{session_id}(:, pcs_to_viz), window);
    
    figure
    
    % generate symmetric axes limits
    ax = axes(); 
    limx = max(abs([sbp25_sm(:, 1); sbp50_sm(:, 1); sbp75_sm(:, 1)]));
    limy = max(abs([sbp25_sm(:, 2); sbp50_sm(:, 2); sbp75_sm(:, 2)]));
    limz = max(abs([sbp25_sm(:, 3); sbp50_sm(:, 3); sbp75_sm(:, 3)]));
    xlim(ax, [-ceil(limx) ceil(limx)]);
    ylim(ax, [-ceil(limy) ceil(limy)]);
    zlim(ax, [-ceil(limz) ceil(limz)]);
    
    view(ax, 3)

    % plotting
    % for color-coding functionality, see
    % https://www.mathworks.com/help/matlab/ref/patch.html#bur94a4-1
    sbps = {sbp25_sm sbp50_sm sbp75_sm};
    if any(contains(ports_to_viz, "p25"))
        x = sbps{1}(:, 1); 
        y = sbps{1}(:, 2);
        z = sbps{1}(:, 3); 
        z(end) = NaN; 
        c = (1 : size(sbps{1}, 1))'; 
        patch(x, y, z, c, 'EdgeColor','interp','LineWidth', 2, 'LineJoin','round');
    end
    if any(contains(ports_to_viz, "p50"))
        x = sbps{2}(:, 1); 
        y = sbps{2}(:, 2);
        z = sbps{2}(:, 3); 
        z(end) = NaN; 
        c = (1 : size(sbps{2}, 1))'; 
        patch(x, y, z, c, 'EdgeColor','interp','LineWidth', 2, 'LineJoin','round');
    end
    if any(contains(ports_to_viz, "p75"))
        x = sbps{3}(:, 1); 
        y = sbps{3}(:, 2);
        z = sbps{3}(:, 3); 
        z(end) = NaN; 
        c = (1 : size(sbps{3}, 1))'; 
        patch(x, y, z, c, 'EdgeColor','interp','LineWidth', 2, 'LineJoin','round');
    end

    % axes labeling and titling
    xlabel("PC" + sprintf('%d', pcs_to_viz(1))); 
    ylabel("PC" + sprintf('%d', pcs_to_viz(2))); 
    zlabel("PC" + sprintf('%d', pcs_to_viz(3))); 
    title("session " + sprintf('%.1s', scores_by_port.id{session_id})); 
    cb = colorbar; 
    cb.Label.String = 'Trial Number';
end