function accept = frame_vote_decision(b, alpha)
% b: 0/1 decisions over L symbols
need = ceil(alpha * numel(b));
accept = (sum(b) >= need);
end
