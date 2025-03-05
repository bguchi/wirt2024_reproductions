% utility function to identify port contingencies based on R/NR history
% must check whether returned `ids` contain repeated elements
function ids = port_ids(Event_timestamps, movmean_window)
    ids = strings(3, 1); 
    hyp = stats(); 

    % filter for p1/p2 events based off associated r/nr codes
    p1_rnr = ismember(Event_timestamps(:, 3), [4 7]);
    p2_rnr = ismember(Event_timestamps(:, 3), [5 8]); 

    % calculate running averages of R/NR history
    mvm1 = movmean(Event_timestamps(p1_rnr, 3), movmean_window); 
    mvm2 = movmean(Event_timestamps(p2_rnr, 3), movmean_window);

    % classify initial port contingency based off how running averages "sandwich" around expected R/NR values 
    % assumes ESP occurs after trial 20
    if hyp(1, 1, 1) - .25*hyp(1, 2, 1) < mean(mvm1(1:20)) & mean(mvm1(1:20)) < hyp(1, 1, 1) + .25*hyp(1, 2, 1)
        ids(1) = "25%"; 
    elseif hyp(2, 1, 1) - .25*hyp(2, 2, 1) < mean(mvm1(1:20)) & mean(mvm1(1:20)) < hyp(2, 1, 1) + .25*hyp(2, 2, 1)
        ids(1) = "50%"; 
    else 
        ids(1) = "75%"; 
    end

    if hyp(1, 1, 2) - .25*hyp(1, 2, 2) < mean(mvm2(1:20)) & mean(mvm2(1:20)) < hyp(1, 1, 2) + .25*hyp(1, 2, 2)
        ids(2) = "25%"; 
    elseif hyp(2, 1, 2) - .25*hyp(2, 2, 2) < mean(mvm2(1:20)) & mean(mvm2(1:20)) < hyp(2, 1, 2) + .25*hyp(2, 2, 2)
        ids(2) = "50%"; 
    else 
        ids(2) = "75%"; 
    end
    
    rs = ["25%", "50%", "75%"]'; 
    ids(3) = rs(find(~ismember(rs, ids), 1)); 

    clear hyp p1_rnr p2_rnr mvm1 mvm2 rs
end

% utility statistics functions
function ev = e_v(vs, ps)
    % does not check size(vs) == size(ps)
    ev = sum(vs .* ps); 
end

function sd = s_d(vs, ps)
    % does not check size(vs) == size(ps)
    ev = e_v(vs, ps); 
    sq_diffs = (vs - ev).^ 2;

    sd = sqrt(e_v(sq_diffs, ps));
    clear ev sq_diffs
end

% generates summary statistics 
% over combinations of reward encodings (4/7, 5/8, 6/9) and contingencies (25/75, 50/50, 75/25)
function hyp_stats = stats()
    hyp_stats = zeros(3, 2, 3); 
    vs = zeros(2, 1, 3);
    p25 = [.25 .75]'; 
    p50 = [.5 .5]'; 
    p75 = [.75 .25]'; 

    for i = 1:3
        vs(:, :, i) = [i + 3, i + 6]';
        hyp_stats(:, :, i) = [e_v(vs(:, :, i), p25), s_d(vs(:, :, i), p25); 
                              e_v(vs(:, :, i), p50), s_d(vs(:, :, i), p50); 
                              e_v(vs(:, :, i), p75), s_d(vs(:, :, i), p75)]; 
    end
    clear vs p25 p50 p75 i
end