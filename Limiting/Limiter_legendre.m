classdef Limiter_legendre < Limiter
    % Base class for all modal limiters, i.e. those which limit Legendre
    % coefficients of the solution polynomials (e.g. hierarchical moment
    % limiters, WENO limiters).
    % In charge of pre- and post-processing DoFs of any basis such that a 
    % derived class implementation will have padded Legendre coefficients 
    % of each troubled element readily available. Will only replace the
    % DOFs of the solution if it detects a change between limited and
    % unlimited versions of its Legendre coefficients.
    %
    properties (Access = protected)
        % For every troubled element, four 3D arrays of Legendre 
        % coefficients are stored in memory (simultaneously).
        % At the time that applyNoSync is called (by default), arrangement
        % is: element (row), coefficient (column), equation (page). 
        coefs,coefs0,coefsL,coefsR
        % Three scalars keep track of the padded 3D array sizes (equation,
        % coefficient and element counts).
        I,J,K
        % The mapping from local characteristic to conservative variables
        % (i.e. right eigenvector matrix) of each troubled element is also
        % kept in memory.
        R
    end
    methods (Abstract, Access = protected)
        % Apply this limiter equation-wise, on local characteristic fields.
        applyNoSync(this)
    end
    methods
        %% Constructor
        function this = Limiter_legendre(varargin)
            % Superclass constructor:
            this = this@Limiter(varargin{:});
        end
        %% Apply (extension)
        function applyStage(this,mesh,solver)
            % Default limiting:
            applyStage@Limiter(this,mesh,solver);
            % Retrieve troubled elements:
            elements = mesh.elements([mesh.elements.isTroubled]);
            % Initialize all the stuff:
            this.initialize(elements,mesh.maxBasisCount)
            % Limit each characteristic component (independently):
            this.applyNoSync;
            % Rewrite limited coefficients into troubled elements:
            this.rewrite(elements);
        end
    end
    methods (Access = protected)
        %% Initialize auxiliary properties
        function initialize(this,elements,N)
            % Intializes various auxiliary variables necessary for any
            % limiting procedure. This approach incurs in a significant 
            % memory consumption penalty, but simplifies the algorithm, 
            % and reduces computational cost overall.
            %
            % Arguments
            %  elements: a row array of elements to be limited.
            %  N: number of Legendre coefficients to compute per element.
            %
            this.I = this.physics.equationCount;
            this.J = N - 1; % second-highest Legendre coefficient to limit (J limits J+1)
            this.K = length(elements);
            % Preallocate arrays of Legendre coefficients and differences:
            this.coefs = zeros(this.I,N,this.K);
            this.coefsL = this.coefs;
            this.coefsR = this.coefs;
            % Preallocate inverse characteristic projection operators:
            this.R = zeros(this.I,this.I,this.K);
            % Populate Legendre coefficients:
            for k = 1:this.K
                this.coefs(:,1:elements(k).dofCount,k) = elements(k).getLegendre;
            end
            % Save unlimited Legendre coefs. in cons. vars. for later:
            this.coefs0 = this.coefs;
            % Populate remaining stuff:
            for k = 1:this.K
                % Get local characteristic projection operators:
                [~,L,this.R(:,:,k)] = this.physics.getEigensystemAt(this.coefs(:,1,k));
                % Apply to Legendre coefficients:
                this.coefs(:,:,k) = L*this.coefs(:,:,k);
                this.coefsL(:,1:elements(k).elementL.dofCount,k) = L*elements(k).elementL.getLegendre;
                this.coefsR(:,1:elements(k).elementR.dofCount,k) = L*elements(k).elementR.getLegendre;
            end
            % Permute for convenience:
            this.coefs = permute(this.coefs,[3,2,1]);
            this.coefsL = permute(this.coefsL,[3,2,1]);
            this.coefsR = permute(this.coefsR,[3,2,1]);
        end
        %% Rewrite element solutions
        function rewrite(this,elements)
            % Replaces the state matrix of an array of elements with the
            % projection to each's basis of the limited Legendre
            % coefficients living in this limiter (projected back to
            % conservative variables).
            %
            % Return limited coefficients to their original arrangement:
            this.coefs = permute(this.coefs,[3,2,1]);
            % Loop over troubled elements:
            for k = 1:this.K
                % Project back to conservative variables:
                this.coefs(:,:,k) = this.R(:,:,k)*this.coefs(:,:,k);
                % Determine which conservative variables have been affected
                % by the limiting (which was done in characteristic ones):
                elements(k).isLimited = abs(this.coefs0(:,:,k) - this.coefs(:,:,k)) > 1e-10;
                % Update all modes of each affected conservative variables:
                i = any(elements(k).isLimited,2);
                elements(k).setLegendre(this.coefs(i,1:elements(k).dofCount,k),i);
            end
        end
    end
end