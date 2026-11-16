% Copyright (c) 2026 Caiming Tan
% SPDX-License-Identifier: MIT
function [selected_a,selected_b,alarm] = ...
    selectCandidateWithMinimalInformation(i,j,a_candidates,b_candidates, ...
    sign_a_i,b_j1,b_j2,criterion,threshold)
% Select an MIAD candidate using one additional component and, when needed,
% the sign of a_i.

if size(a_candidates,1) ~= 3 || size(b_candidates,1) ~= 3 ...
        || size(a_candidates,2) ~= size(b_candidates,2)
    error('MIAD:InvalidCandidateSet', ...
        'Candidate vectors must be 3-by-N arrays with matching N.');
end
if isempty(a_candidates)
    error('MIAD:NoCandidates', 'No feasible attitude candidates were supplied.');
end

candidate_indices = 1:size(a_candidates,2);
if i == j
    if criterion
        discriminator_axis = cyclicIndex(j,1);
        observed_discriminator = b_j1;
    else
        discriminator_axis = cyclicIndex(j,2);
        observed_discriminator = b_j2;
    end
elseif i == cyclicIndex(j,1)
    if criterion
        discriminator_axis = cyclicIndex(j,1);
        observed_discriminator = b_j1;
    else
        candidate_indices = filterBySign(candidate_indices, ...
            a_candidates(i,:),sign_a_i);
        discriminator_axis = cyclicIndex(j,2);
        observed_discriminator = b_j2;
    end
else
    if criterion
        candidate_indices = filterBySign(candidate_indices, ...
            a_candidates(i,:),sign_a_i);
        discriminator_axis = cyclicIndex(j,1);
        observed_discriminator = b_j1;
    else
        discriminator_axis = cyclicIndex(j,2);
        observed_discriminator = b_j2;
    end
end

if isempty(candidate_indices)
    error('MIAD:NoCandidatesAfterSignFilter', ...
        'The sign constraint removed every attitude candidate.');
end
residuals = abs(b_candidates(discriminator_axis,candidate_indices) ...
    - observed_discriminator);
[sorted_residuals,order] = sort(residuals);
selected_index = candidate_indices(order(1));
selected_a = a_candidates(:,selected_index);
selected_b = b_candidates(:,selected_index);
alarm = numel(sorted_residuals) > 1 ...
    && sorted_residuals(2) - sorted_residuals(1) <= threshold;
end

function candidate_indices = filterBySign(candidate_indices,a_i_candidates,sign_a_i)
if sign_a_i > 0
    candidate_indices = candidate_indices(a_i_candidates(candidate_indices) >= 0);
else
    candidate_indices = candidate_indices(a_i_candidates(candidate_indices) <= 0);
end
end
