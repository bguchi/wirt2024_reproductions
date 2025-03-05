function [filtered_and_avgd_iFRs, np_ports] = filter_and_avg_iFRs(session_id)
    load("3portdata/" + session_id, "Event_timestamps", "iFR", "Tmtx");
    start = Tmtx(1); 
    sr = 1 / (Tmtx(2) - start); 

    np_mask = ismember(Event_timestamps(:, 3), 1:3); 
    np_intervals = Event_timestamps(np_mask, 1:2);
    np_ports = Event_timestamps(np_mask, 3); 
    
    filtered_and_avgd_iFRs = zeros(sum(np_mask), size(iFR, 1));
    for trial = 1 : sum(np_mask)
        start_tmtx = floor((np_intervals(trial, 1) - start) * sr); 
        end_tmtx = floor((np_intervals(trial, 2) - start) * sr); 
        filtered_and_avgd_iFRs(trial, :) = mean(iFR(:, start_tmtx:end_tmtx), 2)'; 
    end
    clearvars -except filtered_and_avgd_iFRs np_ports
end