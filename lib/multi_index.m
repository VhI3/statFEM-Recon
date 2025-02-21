function alpha = multi_index(M, p)
% MULTI_INDEX Generates a multi-index sequence for M-dimensional polynomials.
%
%   alpha = MULTI_INDEX(M, p) computes the multi-index sequence for
%   M-dimensional polynomials up to the given order p.
%
% Inputs:
%   M - Number of dimensions (random variables).
%   p - Polynomial order.
%
% Outputs:
%   alpha - Multi-index matrix, where each row represents the indices
%           of a term in the polynomial basis.
%
% Notes:
%   - The output includes all polynomial terms up to order p.
%   - Optimized for simplicity and performance.
%
% Project: statFEM-Recon
% Author: Vahab Narouie
% License: GNU GPL v3.0 (see LICENSE file for details)
% Adapted and optimized from the original by Felipe Uribe, Oct/2014

%% Handle Special Case for M = 1
if M == 1
    alpha = (0:p)'; % Directly return for 1D case
    return;
end

%% Initialize Multi-Index
alpha = zeros(0, M); % Preallocate empty matrix for results

%% Compute Multi-Index Sequence for M > 1
for q = 0:p
    % Generate all combinations of indices that sum to q
    combinations = nchoosek(1:(M + q - 1), M - 1);
    differences = diff([zeros(size(combinations, 1), 1), combinations, (M + q) * ones(size(combinations, 1), 1)], 1, 2) - 1;
    alpha = [alpha; differences]; %#ok<AGROW> % Append new indices
end

end
