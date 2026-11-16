% Copyright (c) 2026 Caiming Tan
% SPDX-License-Identifier: MIT
function [selected_a,selected_b] = selectCandidateWithFullInformation( ...
    i,j,a_candidates,b_candidates,a_i,b_j1,b_j2)
% Select the candidate minimizing residuals in all available components.

candidate_a_i = a_candidates(i,:);
candidate_b_j1 = b_candidates(cyclicIndex(j,1),:);
candidate_b_j2 = b_candidates(cyclicIndex(j,2),:);
cost = (candidate_a_i - a_i).^2 ...
    + (candidate_b_j1 - b_j1).^2 ...
    + (candidate_b_j2 - b_j2).^2;
[~,best_candidate] = min(cost);
selected_a = a_candidates(:,best_candidate);
selected_b = b_candidates(:,best_candidate);
end
