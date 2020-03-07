classdef BDF < Limiter
    % Hierarchical moment limiter by Biswas, Devine and Flaherty (1994).
    % Free of user-defined parameters. Does not guarantee any TVD or TVB 
    % property (not even in the means).
    %
    properties (SetAccess = immutable)
        % Biswas et al., 1994 suggest carrying out 2 limiter sweeps each
        % time the limiter is called. This ensures that the modes used to
        % determine if limiting is required have themselves been limited.
        sweepCount
    end
    properties (Access = protected)
        % All the states of every troubled element and its neighbors are
        % duplicated and stored in the form of Legendre coefficients and
        % left/right-sided differences of characteristic variables. These
        % are 3D arrays in which each row is an element, each column a
        % Legendre coefficient, and each page a characteristic variable.
        % All three have the same number of rows and pages, but coefs has 
        % one more column than the other two.
        coefs,difsL,difsR
        % An additional 3D array is used to store the limited status of
        % every element, component and equation.
        extra
        % The mapping from local characteristic to conservative variables
        % (i.e. right eigenvector matrix) of each troubled element is also
        % stored in memory.
        R
        % Three scalars indicate the (initial) length of the previous:
        I,J,K
    end
    methods
        %% Constructor
        function this = BDF(varargin)
            % Superclass constructor:
            this = this@Limiter(varargin{:});
            % Initialize an input parser:
            p = inputParser;
            p.KeepUnmatched = true;
            addParameter(p,'sweeps',2);
            % Parse the inputs:
            parse(p,varargin{:});
            this.sweepCount = p.Results.sweeps;
        end
        %% Apply (extension)
        function apply(this,mesh,~)
            % Default limiting:
            apply@Limiter(this,mesh);
            % Retrieve troubled elements:
            elements = findobj(mesh.elements,'isTroubled',true)';
            % Perform a number of limiter sweeps:
            for s = 1:this.sweepCount
                % Initialize all the stuff:
                this.initialize(elements)
                % Limit each characteristic component (independently):
                this.applyNoSync;
                % Rewrite limited coefficients into troubled elements:
                this.rewrite(elements);
            end
        end
    end
    methods (Access = protected)
        %% Initialize auxiliary properties
        function initialize(this,elements)
            % Intializes various auxiliary variables necessary for the
            % limiting procedure. This approach incurs in a significant 
            % memory consumption penalty, but simplifies the algorithm, 
            % and reduces computational cost overall.
            %
            %  elements: a row array of elements to be limited.
            %
            % Initialize size variables:
            this.I = this.physics.equationCount;
            this.J = max([elements.dofCount]) - 1; % second-highest Legendre coefficient to limit (J limits J+1)
            this.K = length(elements);
            % Compute coefficient difference weights:
            weights = 1./(2*(1:this.J)-1);
            % Preallocate arrays of Legendre coefficients and differences:
            this.coefs = zeros(this.K,this.J+1,this.I);
            this.difsL = zeros(this.K,this.J,this.I);
            this.difsR = this.difsL;
            % Preallocate inverse characteristic projection operators:
            this.R = zeros(this.I,this.I,this.K);
            % Populate them:
            for k = 1:this.K
                % Get Legendre coefficients:
                aux = elements(k).basis.getLegendre(elements(k));
                % Local characteristic projection:
                [~,L,this.R(:,:,k)] = this.physics.getEigensystemAt(aux(:,1));
                aux = L*aux;
                % Legendre coefficients of neighbors:
                auxL = L*elements(k).edgeL.elementL.basis.getLegendre(elements(k).edgeL.elementL);
                auxR = L*elements(k).edgeR.elementR.basis.getLegendre(elements(k).edgeR.elementR);
                % Reshaping + assignment + padding:
                for page = 1:this.I
                    % Troubled element:
                    this.coefs(k,1:size(aux,2),page) = aux(page,:);
                    % Left-wise differences:
                    this.difsL(k,1:size(aux,2)-1,page) = aux(page,1:end-1);
                    this.difsL(k,1:size(auxL,2)-1,page) = this.difsL(k,1:size(auxL,2)-1,page) - auxL(page,1:end-1);
                    this.difsL(k,:,page) = weights.*this.difsL(k,:,page);
                    % Right-wise differences:
                    this.difsR(k,1:size(aux,2)-1,page) = -aux(page,1:end-1);
                    this.difsR(k,1:size(auxR,2)-1,page) = this.difsR(k,1:size(auxR,2)-1,page) + auxR(page,1:end-1);
                    this.difsR(k,:,page) = weights.*this.difsR(k,:,page);
                end
            end
            % Duplicate the unlimited coefficients for later:
            this.extra = this.coefs;
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
                    % Store limited coefs.:
                    this.coefs(k,j+1,i) = this.minmod(this.coefs(k,j+1,i),this.difsL(k,j,i),this.difsR(k,j,i));
                    % Exclude elements for which there was no change:
                    k(k) = abs(this.extra(k,j+1,i) - this.coefs(k,j+1,i)) > 1e-10;
                    %%% Store limiter activation in characteristic variables:
                    this.difsR(:,j,i) = k;
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
            this.extra = permute(this.extra,[3,2,1]);
            this.coefs = permute(this.coefs,[3,2,1]);
            % Loop over troubled elements:
            for k = 1:this.K
                % Project back to conservative variables:
                this.extra(:,:,k) = this.R(:,:,k)*this.extra(:,:,k);
                this.coefs(:,:,k) = this.R(:,:,k)*this.coefs(:,:,k);
                % Determine which conservative variables have been affected
                % by the limiting (which was done in characteristic ones):
                elements(k).isLimited = abs(this.extra(:,:,k) - this.coefs(:,:,k)) > 1e-10;
                % Replace all modes of affected conservative variables:
                i = any(elements(k).isLimited,2);
                elements(k).basis.setLegendre(elements(k),this.coefs(i,:,k),i);
            end
        end
    end
    methods (Static, Access = protected)
        %% Set limited coefficients
        function setLimitedLegendre(element,coefs)
            % Find all rows of the given element's state matrix that have
            % been limited and replace them by the projection to the 
            % element's basis of the given matrix of Legendre coefficients.
            %
            rows = any(element.isLimited,2);
            element.basis.setLegendre(element,coefs(rows,:),rows);
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