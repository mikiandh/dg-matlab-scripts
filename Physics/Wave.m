classdef Wave < Physics
    % Wave equation physics. Hyperbolic system of two linear equations with
    % constant coefficients. Alternative formulation of the scalar 2nd 
    % order wave equation.
    properties (Constant, Hidden)
        equationCount = 2
        controlVars = 1:2
    end
    properties (SetAccess = immutable)
        waveSpeed = 1
    end
    properties (Access = protected) %%% might break something, check later
        jacobian % i.e. A matrix
        jacobianL % i.e. A^+ or A^- matrix (depending on sign of waveSpeed)
        jacobianR % idem
        eigenvectors
        invEigenvectors
    end
    methods
        %% Constructor
        function this = Wave(c)
            if nargin
                validateattributes(c,"numeric",{'finite'})
                this.waveSpeed = c;
            end
            this.jacobian = [0 this.waveSpeed^2; 1 0];
            this.jacobianL = 0.5*[abs(this.waveSpeed) this.waveSpeed^2; 1 abs(this.waveSpeed)];
            this.jacobianR = [-1 1; 1 -1].*this.jacobianL;
            this.eigenvectors = [this.waveSpeed -this.waveSpeed; 1 1];
            this.invEigenvectors = [1 this.waveSpeed; -1 this.waveSpeed]/(2*this.waveSpeed);
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
        %% Sample Jacobian matrix
        function A = getJacobianAt(this,~)
            % Returns the Jacobian matrix corresponding to a given state
            % vector.
            %
            A = this.jacobian;
        end
        %% Jacobian eigen-decomposition
        function [D,L,R] = getEigensystemAt(this,~,~)
            % Returns the eigenvalue and eigenvector matrices evaluated at,
            % either:
            %
            % A) the given state vector
            % B) the "generalized Roe average" between two given states
            %
            D = this.waveSpeed*[1 0; 0 -1];
            L = this.invEigenvectors;
            R = this.eigenvectors;
        end
        %% Advection velocity
        function v = getVelocityAt(this,~)
            v = this.waveSpeed;
        end
        %% Information (extension)
        function info = getInfo(this)
            info = sprintf('%s, c = %g',this.getInfo@Physics,this.waveSpeed);
        end
    end
end