function [vectSurf] = makeSurf(elementNodes, vect)
% MAKESURF Generates surface data for visualization.
%
%   [vectSurf] = MAKESURF(elementNodes, vect) organizes the input vector
%   into a surface format based on the element connectivity.
%
% Inputs:
%   elementNodes - (nElements x 4) Matrix of element connectivity, where each
%                  row lists the node indices for a quadrilateral element.
%   vect         - (nNodes x 1) Vector of values associated with each node.
%
% Outputs:
%   vectSurf     - (4 x nElements) Matrix of values organized for visualization.
%                  Each column corresponds to an element's surface data.
%
% Notes:
%   - The function assumes quadrilateral elements (4 nodes per element).
%   - Designed for use in finite element visualization tasks.
%
% Project: statFEM-Recon
% Author: Vahab Narouie
% License: GNU GPL v3.0 (see LICENSE file for details)

%% Ensure Input Vector is Column-Oriented
if size(vect, 2) ~= 1
    vect = vect'; % Transpose to column vector if needed
end

%% Map Node Values to Element Surfaces
nElements = size(elementNodes, 1); % Number of elements
nodeIndices = elementNodes(1:nElements, :); % Connectivity indices
vectSurf = vect(nodeIndices', 1); % Map node values to element surfaces
vectSurf = reshape(vectSurf, [4, nElements]); % Reshape for quadrilateral format
end
