classdef BDF < Limiter
    % Hierarchical moment limiter by Biswas, Devine and Flaherty (1994).
    % Free of user-defined parameters. Does not guarantee any TVD or TVB 
    % property (not even in the means).
    %
    properties (Access = protected)
        % For every troubled element, four 3D arrays are stored in memory
        % (simultaneously).
        coefs,coefs_backup,difsL,difsR
        % Three scalars keep track of the aforementioned 3D array sizes.
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
        function apply(this,mesh,~)
            % Default limiting:
            apply@Limiter(this,mesh);
            % Retrieve troubled elements:
            elements = findobj(mesh.elements,'isTroubled',true)';
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
            % Compute coefficient difference weights:
            weights = 1./(2*(1:this.J)-1);
            % Preallocate arrays of Legendre coefficients and differences:
            this.coefs = zeros(this.I,N,this.K);
            this.difsL = zeros(this.I,this.J,this.K);
            this.difsR = this.difsL;
            % Preallocate inverse characteristic projection operators:
            this.R = zeros(this.I,this.I,this.K);
            % Populate Legendre coefficients:
            for k = 1:1:this.K
                this.coefs(:,1:elements(k).dofCount,k) = elements(k).basis.getLegendre(elements(k));
            end
            % Duplicate original Legendre coefficients "just in case":
            this.coefs_backup = this.coefs;
            % Populate the remaining variables:
            for k = 1:this.K
                % Get local characteristic projection operators:
                [~,L,this.R(:,:,k)] = this.physics.getEigensystemAt(this.coefs(:,1,k));
                % Apply local characteristic projection operator:
                this.coefs(:,:,k) = L*this.coefs(:,:,k);
                % Compute left/right-sided finite differences (zero-padded):
                j = 1:elements(k).edgeL.elementL.dofCount-1;
                this.difsL(:,j,k) = L*elements(k).edgeL.elementL.basis.getLegendre(elements(k).edgeL.elementL,j);
                this.difsL(:,:,k) = weights.*(this.coefs(:,1:end-1,k) - this.difsL(:,:,k));
                j = 1:elements(k).edgeR.elementR.dofCount-1;
                this.difsR(:,j,k) = L*elements(k).edgeR.elementR.basis.getLegendre(elements(k).edgeR.elementR,j);
                this.difsR(:,:,k) = weights.*(this.difsR(:,:,k) - this.coefs(:,1:end-1,k));
            end
            % Permute them for convenience:
            this.coefs = permute(this.coefs,[3,2,1]);
            this.difsL = permute(this.difsL,[3,2,1]);
            this.difsR = permute(this.difsR,[3,2,1]);
        end
        %% Apply equation-wise
        function applyNoSync(this)
            % Applies the minmod limiter asynchronously (i.e. on each
            % characteristic variable separately).
            %
            % Loop over characteristic variables:
            for i = 1:this.I
                k = true(this.K,1); % flag every element for limiting
                % Limit hierarchically (top to bottom):
                for j = this.J:-1:1
                    % Get limited coefs.:
                    aux = this.coefs(k,j+1,i);
                    this.coefs(k,j+1,i) = this.minmod(this.coefs(k,j+1,i),this.difsL(k,j,i),this.difsR(k,j,i));
                    % Update list of elements flagged for limiting:
                    k(k) = abs(this.coefs(k,j+1,i) - aux) > 1e-10 | ~aux;
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
            % Permute coefficients to the conventional arrangement:
            this.coefs = permute(this.coefs,[3,2,1]);
            % Loop over troubled elements:
            for k = 1:this.K
                % Project back to conservative variables:
                this.coefs(:,:,k) = this.R(:,:,k)*this.coefs(:,:,k);
                % Determine which conservative variables have been affected
                % by the limiting (which was done in characteristic ones):
                elements(k).isLimited = abs(this.coefs_backup(:,:,k) - this.coefs(:,:,k)) > 1e-10;
                % Replace all modes of affected conservative variables:
                i = any(elements(k).isLimited,2);
                elements(k).basis.setLegendre(elements(k),this.coefs(i,1:elements(k).dofCount,k),i);
            end
        end
    end
    methods (Static)
        %% Information
        function info = getInfo
            info = 'BDF';
        end
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