function BVP = tensionBar_1D_P_matrix(BVP)
% TENSIONBAR_1D_P_MATRIX
% Constructs the projection matrix P for a 1D tension bar.
%
% Inputs:
%   BVP - Boundary Value Problem structure containing mesh and observation data.
%
% Outputs:
%   P - Projection matrix mapping FEM nodes to sensor locations.
%
% Project: statFEM-Recon
% Author: Vahab Narouie
% License: GNU GPL v3.0 (see LICENSE file for details)

%% Step 1: Extract Data from BVP
senCoor = BVP.obs.senCoor; % Sensor positions
nElm = BVP.mesh.nElm; % Number of elements
nodeCoordinates = BVP.mesh.nodeCoordinates;
numberNodes = BVP.mesh.numberNodes;
elementNodes = BVP.fem.elementNodes;
nsen = BVP.obs.nsen;
% Element-wise coordinates
elementCoordinates = nodeCoordinates(elementNodes);

% Initialize projection matrix
P = zeros(nsen, numberNodes);

%% Step 2: Compute Projection Matrix
for i = 1:length(senCoor)
    
    for j = 1:nElm
        % Check if senCoor point is within the element
        if (senCoor(i) >= elementCoordinates(j, 1)) && (senCoor(i) <= elementCoordinates(j, 2))
            % Map to local coordinate ξ in [-1,1]
            xi = (2 * senCoor(i) - sum(elementCoordinates(j, :))) / ...
                (elementCoordinates(j, 2) - elementCoordinates(j, 1));
            
            % Compute shape functions
            N = Lagrange1D(xi, 1, [-1, 1]);
            
            % Assign values to projection matrix
            P(i, elementNodes(j, :)) = N;
            break; % Found the element, no need to check further
        end
        
    end
    
end

%% Assign Back to BVP
BVP.obs.P = P;
BVP.obs.P_active = P(1:nsen, 2:end);
disp('6. Projection Matrix for 1D Tension Bar Completed.');
end
