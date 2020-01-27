classdef Wave < Physics
    properties
        waveSpeed
        jacobian % i.e. A matrix
        jacobianL % i.e. A^+ or A^- matrix (depending on sign of waveSpeed)
        jacobianR % idem
        eigenvectors
        invEigenvectors
    end
    methods
        %% Constructor
        function wave = Wave(c)
            if nargin == 0
                c = 1;
            end
            wave.equationCount = 2;
            wave.waveSpeed = c;
            wave.jacobian = [0 c^2; 1 0];
            wave.jacobianL = 0.5*[abs(c) c^2; 1 abs(c)];
            wave.jacobianR = 0.5*[-abs(c) c^2; 1 -abs(c)];
            wave.eigenvectors = [c -c; 1 1];
            wave.invEigenvectors = [1 c; -1 c]/(2*c);
        end
        %% Convection operator
        function convection = getConvectionMatrix(this,~,basis)
            % Returns the discrete convection operator to be applied to the
            % modal coefficients of the discretizated solution.
            %
            % Arguments
            %  states: (unused) row array of basis function coefficients
            %  basis: finite-dimensional discretization basis
            % Output
            %  convection: discrete convection operator
            error('Vector physics not supported.')
        end
        %% Flux function
        function flux = flux(this,states)
            flux = this.jacobian*states;
        end
        %% Riemann solver (exact)
        function [flux,S] = riemannFlux(this,stateL,stateR)
            flux = this.jacobianL*stateL + this.jacobianR*stateR;
            S = this.waveSpeed*[1 -1]';
        end
        %% Outlet boundary conditions
%         function applyBoundaryConditionsWeakly(this,mesh)
%             externalState = zeros(2,1);
%             % Right-most boundary:
%             mesh.elements(end).riemannR =...
%                 this.riemannFlux(mesh.elements(end).stateR,externalState);
%             mesh.edges{end}.computeTimeDeltas(-this.waveSpeed);
%             % Left-most boundary:
%             mesh.elements(1).riemannL =...
%                 -this.riemannFlux(externalState,mesh.elements(1).stateL);
%             mesh.edges{1}.computeTimeDeltas(this.waveSpeed);
%         end
        %% Reflecting boundary conditions
        function applyBoundaryConditions(this,mesh)
            % Right-most boundary:
            externalState = [mesh.elements(end).stateR(1); -mesh.elements(end).stateR(2)];
            mesh.elements(end).riemannR = this.riemannFlux(mesh.elements(end).stateR,externalState);
            mesh.edges{end}.computeTimeDeltas(-this.waveSpeed);
            % Left-most boundary:
            externalState = [mesh.elements(1).stateL(1); -mesh.elements(1).stateL(2)];
            mesh.elements(1).riemannL =...
                -this.riemannFlux(externalState,mesh.elements(1).stateL);
            mesh.edges{1}.computeTimeDeltas(this.waveSpeed);
        end
        %% Convert state vector(s) to characteristic variable vector(s)
        function states = stateToEigenstate(this,states)
            states = this.invEigenvectors*states;
        end
        %% Convert characteristic variable vector(s) to state vector(s)
        function states = eigenstateToState(this,states)
            states = this.eigenvectors*states;
        end
        %% Return 3D arrays of left and right eigenvectors
        function [L,R] = getEigenvectors(this,meanStates)
            % Returns the eigenvector matrices corresponding to each
            % mean state in the input, including left-averaged and
            % right-averaged (left) eigenvectors.
            %
            % meanStates(i,j,k)
            % L.{'left','none','right'}(i,j,k)
            % R.{'left','none','right'}(i,j,k)
            %  struct fields: left-averaged, not averaged, right-averaged
            %  i: system dimension
            %  j: eigenvector index
            %  k: element
            K = size(meanStates,3) - 2; % 1st and last states are excluded
            L.left = repmat(this.invEigenvectors,1,1,K); % Jacobian matrix is constant in the wave equation
            L.none = repmat(this.invEigenvectors,1,1,K); % Jacobian matrix is constant in the wave equation
            L.right = repmat(this.invEigenvectors,1,1,K); % Jacobian matrix is constant in the wave equation
            R.left = repmat(this.eigenvectors,1,1,K); % Jacobian matrix is constant in the wave equation
            R.none = repmat(this.eigenvectors,1,1,K); % Jacobian matrix is constant in the wave equation
            R.right = repmat(this.eigenvectors,1,1,K); % Jacobian matrix is constant in the wave equation
        end
        %% Return left and right eigenvectormatrices for given mean states
        function [L,R] = getEigenvectorsAt(this,meanStates)
            % Returns the eigenvector matrices evaluated at the mean state
            % passed as input.
            %
            % meanStates(i,j,k)
            % L(i,j,k)
            % R(i,j,k)
            %  i: system dimension
            %  j: eigenvector index
            %  k: element
            %
            K = size(meanStates,3);
            L = repmat(this.invEigenvectors,1,1,K);
            R = repmat(this.eigenvectors,1,1,K);
        end
    end
    methods (Static)
        %% Check for positivity in certain quantities
        function PASS = checkWithinPhysicalBounds(~,~)
            PASS = 1;
        end
    end
end