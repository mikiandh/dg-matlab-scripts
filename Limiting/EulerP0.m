classdef EulerP0 < Limiter_legendre
    % Failsafe limiter for Euler equations (Godunov scheme version).
    %
    % Sets linear and above Legendre coefficients (in characteristic
    % variables) to zero, on those elements in which it detects unphysical
    % states.
    %
    methods
        %% Constructor
        function this = EulerP0(varargin)
            % Superclass constructor:
            this = this@Limiter_legendre(varargin{:});
        end
        %% Apply (override)
        function applyStage(this,mesh,solver)
            % Default limiting:
            this.applyStage@Limiter(mesh,solver);
            % Get elements with invalid states:
            elements = this.getFailedElements(mesh);
            if ~isempty(elements)
                % Keep only those which are also troubled:
                isTroubled = [elements.isTroubled];
                elements(~isTroubled(:,:,this.priority)) = [];
                % Extract their Legendre coefficients and all that stuff:
                this.initialize(elements,mesh.maxBasisCount)
                % Limit each characteristic component (independently):
                this.applyNoSync
                % Rewrite limited coefficients:
                this.rewrite(elements)
            end
        end
        %% Apply (extension)
        function applyInitial(this,mesh,solver)
            validateattributes(solver.physics,{'Euler'},{})
            this.applyInitial@Limiter_legendre(mesh,solver);
        end
    end
    methods (Access = protected)
        %% Apply equation-wise (override)
        function applyNoSync(this)
            % Last-resort: keep only the mean value.
            this.coefs(:,2:end,:) = 0;
        end
        %% Unphysical states (interior and edges, recommended)
        function elements = getFailedElements(this,mesh)
            % Extrapolate state at all edges:
            mesh.elements.interpolateStateAtEdges
            % Search for troublesome state pairs at edges:
            isTroubledEdge = true(1,mesh.edgeCount);
            for k = 1:mesh.edgeCount
                try
                    states = [mesh.edges(k).elementL.stateR mesh.edges(k).elementR.stateL];
                    vars = Euler.allPrimitiveVarsFromState(states);
                catch
                    continue % if this fails, the edge is troubled for sure
                end
                if 5*sum(vars(5,:)) <= diff(vars(2,:))
                    continue % unsafe, vacuum generation (Toro, eq. 4.40)
                end
                vars(2,:) = [];
                if any(vars(:) < 0)
                    continue % invalid negative stuff
                end
                % We got all the way down here, so no troubles were found:
                isTroubledEdge(k) = false;
            end
            % Search for invalid states within troubled elements:
            isTroubledElement = [mesh.elements.isTroubled];
            isTroubledElement = isTroubledElement(:,:,this.priority);
            for k = find(isTroubledElement)
                if isTroubledEdge(k) || isTroubledEdge(k+1)
                    % An element can only be free of troubles if both its
                    % edges are so too.
                    continue
                end
                coords = mesh.elements(k).getGaussCoords;
                states = mesh.elements(k).interpolateStateAtCoords(coords');
                try
                    vars = Euler.allPrimitiveVarsFromState(states);
                catch
                    continue % if this fails, element is troubled
                end
                % Ignore velocity:
                vars(2,:) = [];
                % If any of the remaining is negative, we are in trouble:
                isTroubledElement(k) = any(vars(:) < 0);
            end
            % Return those elements with invalid states:
            elements = mesh.elements(isTroubledElement);
        end
        %% Unphysical states (only checks edges)
        function elements = getFailedElements_edges(this,mesh)
            % Extrapolate state at all edges:
            mesh.elements.interpolateStateAtEdges
            % Search for troublesome state pairs at edges:
            isTroubledEdge = true(1,mesh.edgeCount);
            for k = 1:mesh.edgeCount
                try
                    states = [mesh.edges(k).elementL.stateR mesh.edges(k).elementR.stateL];
                    vars = Euler.allPrimitiveVarsFromState(states);
                catch
                    continue % if this fails, the edge is troubled for sure
                end
                if 5*sum(vars(5,:)) <= diff(vars(2,:))
                    continue % unsafe, vacuum generation (Toro, eq. 4.40)
                end
                vars(2,:) = [];
                if any(vars(:) < 0)
                    continue % invalid negative stuff
                end
                % We got all the way down here, so no troubles were found:
                isTroubledEdge(k) = false;
            end
            % Search for invalid states within troubled elements:
            isTroubledElement = [mesh.elements.isTroubled];
            isTroubledElement = isTroubledElement(:,:,this.priority);
            k = find(isTroubledElement(:,:,this.priority));
            isTroubledElement(k) = isTroubledEdge(k) || isTroubledEdge(k+1);
            % Return those elements with invalid states:
            elements = mesh.elements(isTroubledElement);
        end
    end
end