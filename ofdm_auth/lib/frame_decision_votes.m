function accept = frame_decision_votes(votes, L, alpha)
% votes: integer number of 1's over L symbols
accept = (votes >= ceil(alpha * L));
end
