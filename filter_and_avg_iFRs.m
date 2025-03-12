function [filtered_and_avgd_iFRs, np_ports] = filter_and_avg_iFRs(session_id)
    % FILTER_AND_AVG_IFRS Return processed iFRs and port entry sequence for a given session. 
    %    Args:
    %        - session_id: (int) index for corresponding `3portdata` file
    %    Out: 
    %        - filtered_and_avgd_iFRs: (num_trials x num_cells) iFRs, as stored in 
    %          corresponding `3portdata` file, filtered and averaged over NP durations
    %        - np_ports: (num_trials x 1) column vector of port entered (1, 2, or 3) 
    %          for trial i
    load("3portdata/" + session_id, "Event_timestamps", "iFR", "Tmtx");
    start = Tmtx(1); 
    sr = 1 / (Tmtx(2) - start); 

    np_mask = ismember(Event_timestamps(:, 3), 1:3); 
    np_intervals = Event_timestamps(np_mask, 1:2);
    np_ports = Event_timestamps(np_mask, 3); 
    
    filtered_and_avgd_iFRs = zeros(sum(np_mask), size(iFR, 1));
    % iFR uses Tmtx-corresonding indices; elements of Tmtx spaced by sampling rate (sr) 
    for trial = 1 : sum(np_mask)
        start_tmtx = floor((np_intervals(trial, 1) - start) * sr); 
        end_tmtx = floor((np_intervals(trial, 2) - start) * sr); 
        filtered_and_avgd_iFRs(trial, :) = mean(iFR(:, start_tmtx:end_tmtx), 2)'; 
    end
end