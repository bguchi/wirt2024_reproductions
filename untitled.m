listing = dir("3portdata"); 
listing = struct2table(listing);
session_ids = listing.name(~listing.isdir)

good_sessions = readmatrix("good_sessions.csv"); 
session_ids = string(session_ids(good_sessions));
clear listing good_sessions
%%