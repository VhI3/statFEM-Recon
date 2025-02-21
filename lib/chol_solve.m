function x = chol_solve(L, b)
% CHOL_SOLVE Solves a system of linear equations using Cholesky decomposition.
%   x = CHOL_SOLVE(L, b) solves the system of equations A*x = b
%   where A = L*L' and L is a lower triangular matrix obtained
%   from the Cholesky decomposition of A.
%
% Inputs:
%   L - Lower triangular matrix from Cholesky decomposition (must satisfy istril(L)).
%   b - Right-hand side vector or matrix.
%
% Outputs:
%   x - Solution to the system of equations A*x = b.
%
% Example usage:
%   L = chol(A, 'lower'); % Compute lower triangular Cholesky factor
%   x = chol_solve(L, b); % Solve the system
%
% Project: statFEM-Recon
% Author: Vahab Narouie
% License: GNU GPL v3.0 (see LICENSE file for details)

% Check if L is a lower triangular matrix
if istril(L)
  % Forward substitution: Solve L*y = b
  y = L \ b;
  
  % Backward substitution: Solve L'*x = y
  x = L' \ y;
else
  % Raise an error if L is not lower triangular
  error('L must be a lower triangular matrix.');
end

end

% function x = chol_solve(L,b)
% if istril( L )
%     y = L\b;
%     x = L'\y;
% else
%     msg = 'L has to be lower trangular.';
%     error(msg)
% end
% end
