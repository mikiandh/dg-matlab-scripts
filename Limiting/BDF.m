classdef BDF < Limiter
    % Hierarchical moment limiter by Biswas, Devine and Flaherty (1994).
    % Free of user-defined parameters. Does not guarantee any TVD or TVB 
    % property (not even in the means).
    %
    properties (Access = protected)
        % For every troubled element, four 3D arrays of Legendre 
        % coefficeints are stored in memory (simultaneously).
        coefs,coefs0,coefsL,coefsR
        % Three scalars keep track of the padded 3D array sizes.
        I,J,K
        % The mapping from local characteristic to conservative variables
        % (i.e. right eigenvector matrix) of each troubled element is also
        % kept in memory.
        R
    end
    methods
        %% Constructor
        function this = BDF(varargin)
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
            % Intializes various auxiliary variables necessary for the
            % limiting procedure. This approach incurs in a significant 
            % memory consumption penalty, but simplifies the algorithm, 
            % and reduces computational cost overall.
            %
            % Arguments
            %  elements: a row array of elements to be limited.
            %  N: number of Legendre coefficients to compute per element.
            %
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
                this.coefs(:,1:elements(k).dofCount,k) = elements(k).basis.getLegendre(elements(k));
            end
            % Save unlimited Legendre coefs. in cons. vars. for later:
            this.coefs0 = this.coefs;
            % Populate remaining stuff:
            for k = 1:this.K
                % Get local characteristic projection operators:
                [~,L,this.R(:,:,k)] = this.physics.getEigensystemAt(this.coefs(:,1,k));
                % Apply to Legendre coefficients:
                this.coefs(:,:,k) = L*this.coefs(:,:,k);
                this.coefsL(:,1:elements(k).elementL.dofCount,k) = L*elements(k).elementL.basis.getLegendre(elements(k).elementL);
                this.coefsR(:,1:elements(k).elementR.dofCount,k) = L*elements(k).elementR.basis.getLegendre(elements(k).elementR);
            end
            % Permute for convenience:
            this.coefs = permute(this.coefs,[3,2,1]);
            this.coefsL = permute(this.coefsL,[3,2,1]);
            this.coefsR = permute(this.coefsR,[3,2,1]);
        end
        %% Apply equation-wise
        function applyNoSync(this)
            % Applies the minmod limiter asynchronously (i.e. on each
            % characteristic variable separately).
            %
            % Loop over characteristic variables:
            for i = 1:this.I
                % Flag every element for limiting:
                k = true(this.K,1);
                % Limit hierarchically (top to bottom):
                for j = this.J:-1:1
                    % Get limited coefs.:
                    coefs_min = this.coefs(k,j+1,i);
                    this.coefs(k,j+1,i) = this.minmod(...
                        (2*j-1).*this.coefs(k,j+1,i),...
                        this.coefs(k,j,i) - this.coefsL(k,j,i),...
                        this.coefsR(k,j,i) - this.coefs(k,j,i))./(2*j-1);
                    % Update list of elements flagged for limiting:
                    k(k) = abs(this.coefs(k,j+1,i) - coefs_min) > 1e-10 | ~coefs_min;
                end
            end
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
                elements(k).basis.setLegendre(elements(k),this.coefs(i,1:elements(k).dofCount,k),i);
            end
        end
    end
    methods (Static)
        %% Minmod operator (Osher, 1984)
        function A = minmod(A,B,C)
            % Minmod function for 3 matrix inputs. Operates entry-wise.
            % "Inlined" for speed. All three matrices are assumed to have 
            % consistent sizes.
            %
            A = (sign(A) == sign(B) & sign(A) == sign(C)).*sign(A).*min(abs(A),min(abs(B),abs(C)));
        end
    end
end