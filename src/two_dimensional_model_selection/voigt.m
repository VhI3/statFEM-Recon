function V = voigt(M)
% VOIGT Transform a symmetric second-order matrix into Voigt notation.
%
% This function converts a 2D (2x2) or 3D (3x3) symmetric tensor into
% its Voigt notation representation, a compact column vector.
%
% Syntax:
%   V = voigt(M)
%
% Input:
%   M - Symmetric matrix (2x2 or 3x3).
%
% Output:
%   V - Voigt column vector:
%       - For 2D (2x2): [M11; M22; M12].
%       - For 3D (3x3): [M11; M22; M33; M12; M23; M13].
%
% Project: statFEM-Recon
% Author: Vahab Narouie, TU-Braunschweig, 2024
% License: GNU GPL v3.0 (see LICENSE file for details)

% Check if M is a square matrix
if ~ismatrix(M) || size(M, 1) ~= size(M, 2)
    error('Input must be a square matrix (2x2 or 3x3).');
end

% Determine size of the input matrix
n = size(M, 1);

if n == 2
    % 2D Voigt notation
    V = [M(1, 1); M(2, 2); M(1, 2)];
elseif n == 3
    % 3D Voigt notation
    V = [M(1, 1); M(2, 2); M(3, 3); M(1, 2); M(2, 3); M(1, 3)];
else
    error('Only 2x2 or 3x3 matrices are supported.');
end

end
