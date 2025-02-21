function [alpha, Psi_s, Psi_p, PsiSqNorm, P] = Hermite_PC(M, p_order)
% HERMITE_PC Computes Hermite Polynomial Chaos basis.
%   [alpha, Psi_s, Psi_p, PsiSqNorm, P] = HERMITE_PC(M, p_order)
%   generates the Polynomial Chaos (PC) basis using Hermite polynomials.
%
% Inputs:
%   M       - Number of random variables (dimensions).
%   p_order - Polynomial order.
%
% Outputs:
%   alpha       - Multi-index for the polynomial basis.
%   Psi_s       - Symbolic representation of the PC basis.
%   Psi_p       - Polynomial coefficients of the PC basis.
%   PsiSqNorm   - Squared norm of each polynomial in the basis.
%   P           - Total number of basis polynomials.
%
% Project: statFEM-Recon
% Author: Vahab Narouie
% License: GNU GPL v3.0 (see LICENSE file for details)
% Originally created by Felipe Uribe (Oct/2014)
% Changes and enhancements: Added comments, cleaned code structure,
%                           updated formatting, and improved readability.

%% Calculate the basis size of Psi
P = 1;

for s = 1:p_order
    P = P + (1 / factorial(s)) * prod(M + (0:s - 1));
end

%% Calculate 1D Hermite polynomials using recurrence relation
% Symbolic representation
syms xi;
He_s = cell(p_order + 1, 1);
He_s{1} = sym(1); % H_0(x) = 1
He_s{2} = xi; % H_1(x) = x

for j = 2:p_order
    He_s{j + 1} = expand(xi * He_s{j} - (j - 1) * He_s{j - 1});
end

% Polynomial coefficients representation
He_p = cell(p_order + 1, 1);
He_p{1} = 1; % H_0(x) = 1
He_p{2} = [1 0]; % H_1(x) = x

for n = 2:p_order
    He_p{n + 1} = [He_p{n} 0] - (n - 1) * [0 0 He_p{n - 1}]; % Recurrence relation
end

%% Define the number of random variables (RVs)
x = cell(1, M);
H_s = cell(p_order + 1, M); % Hermite polynomials for each dimension (symbolic)
H_p = cell(p_order + 1, M); % Hermite polynomials for each dimension (coefficients)

for j = 1:M
    x{j} = sym(sprintf('xi_%d', j)); % Define symbolic variables xi_j
    
    for i = 1:p_order + 1
        H_s{i, j} = subs(He_s{i}, xi, x{j});
        H_p{i, j} = He_p{i};
    end
    
end

%% Compute M-dimensional PC basis
Psi_s = cell(P, 1); % Symbolic version
Psi_p = cell(P, 1); % Polynomial coefficients version
alpha = multi_index(M, p_order); % Create the multi-index for polynomial orders

for i = 2:P + 1
    mult_s = 1; % Initialize symbolic product
    mult_p = 1; % Initialize polynomial product
    
    for j = 1:M
        mult_s = mult_s * H_s{alpha(i - 1, j) + 1, j};
        mult_p = conv(mult_p, H_p{alpha(i - 1, j) + 1, j});
    end
    
    Psi_s{i - 1} = mult_s;
    Psi_p{i - 1} = mult_p;
end

%% Calculate the squared norm of the polynomials
PsiSqNorm = prod(factorial(alpha), 2);

end
