classdef DGIGA_AFC < Bspline
    properties
        lumpedMassMatrixDiagonal % non-zero entries of the lumped mass matrix (in reference element space)
        mode2fluxes % connectivity matrix between a control point (row: control point index) and all other control points within shared support, excluding itself (column: edge linear index)
        diffusionFun = @DGIGA_AFC.diffusionRoeHartenHyman1;
    end
    methods
        %% Constructor
        function this = DGIGA_AFC(varargin)
            this@Bspline(varargin{:});
            if nargin > 0
                this.lumpedMassMatrixDiagonal = full(sum(this.massMatrix,1));
            end
        end
        %% Instantiate from prototype
        function this = clone(prototype,degree)
            this = DGIGA_AFC(prototype.knots,degree,prototype.smoothness);
        end
        %% DGIGA-AFC operator (low-order predictor)
        function computeResiduals(this,element,physics)
            element.computeFluxesFromStates(physics)
            element.residuals = ...
                element.fluxes*this.gradientMatrix...
                - element.riemannR.*this.right'...
                - element.riemannL.*this.left';
            this.diffuseResiduals(element,physics)
            element.residuals = 2/element.dx*element.residuals./this.lumpedMassMatrixDiagonal;
        end
        %% L2 projection (lumped or constrained)
        function project(this,element,fun0)
            % Lumped projection of a function into a finite-dimensional
            % space using a DGIGA basis. Initializes any quantities that
            % the AFC limiter requires to correct the result to its
            % constrained counterpart.
            %
            % Unconstrained L2 projection:
            this.project@Bspline(element,fun0);
            % Store unconstrained projection (high order candidate):
            element.residuals = element.states;
            % Override solution with the lumped projection (predictor):
            element.states = element.states*element.basis.massMatrix./element.basis.lumpedMassMatrixDiagonal;
            % Initialize diffusion matrix blocks to zero:
            [I,N] = size(element.states);
            rj = element.basis.edges(3,:);
            element.diffusions = cell(N);
            element.diffusions(rj) = {sparse(I,I)};
            % Preallocate antidiffusive fluxes:
            element.antidiffusiveFluxes = spalloc(I,N^2,I*length(rj));
        end
    end
    methods (Access = protected)
        %% Low-order predictor residuals
        function diffuseResiduals(this,element,physics)
            % Applies AFC diffusion to the residual matrix, a la Möller &
            % Jaeschke, 2018 (eq. 16). Saves the diffusion block matrix
            % into the given element as a cell array.
            %
            % Precompute the "gradient differences" (Kuzmin et al. 2012, eq. 27):
            eAbs = .5*abs(element.basis.gradientMatrix - element.basis.gradientMatrix');
            % Loop over edges:
            for edge = element.basis.edges
                % Aliases:
                r = edge(1); j = edge(2);
                % Compute a diffusion matrix block:
                element.diffusions{r,j} = eAbs(r,j)*this.diffusionFun(element.states(:,r),element.states(:,j),physics);
                % Accumulate diffusion into the residuals:
                element.residuals(:,r) = element.residuals(:,r) + element.diffusions{r,j}*(element.states(:,j) - element.states(:,r));
            end
        end
    end
    methods (Static)
        %% AFC diffusion a la Kuzmin et al. 2012, eq. 42
        function D = diffusionRoe(state1,state2,physics)
            % Computes a tensorial numerical diffusion matrix block via AFC
            % a la Kuzmin et al. 2012, eq. 42. No entropy fix.
            %
            [D,L,R] = physics.getEigensystemAt(state1,state2);
            D = R*abs(D)*L;
        end
        %% AFC diffusion: Roe + Harten-Hyman(1)
        function D = diffusionRoeHartenHyman1(state1,state2,physics)
            % Computes a tensorial numerical diffusion matrix block via AFC
            % a la Kuzmin et al. 2012, eq. 42; additionally, applies the
            % 1st entropy fix of Harten and Hyman a la Kermani and Plett 
            % 2001, eq. 8.
            %
            [D,L,R] = physics.getEigensystemAt(state1,state2);
            D12 = diag(D);
            D = abs(D12);
            % Fix:
            D1 = diag(physics.getEigensystemAt(state1));
            D2 = diag(physics.getEigensystemAt(state2));
            epsilon = max(0,max(D12-D1,D2-D12));
            i = D < epsilon;
            D(i) = (D(i).^2+epsilon(i).^2)./(2*epsilon(i));
            % AFC:
            D = R*(D.*L);
        end
        %% AFC diffusion: Roe + Harten-Hyman(2)
        function D = diffusionRoeHartenHyman2(state1,state2,physics)
            % Computes a tensorial numerical diffusion matrix block via AFC
            % a la Kuzmin et al. 2012, eq. 42; additionally, applies the
            % 2nd entropy fix of Harten and Hyman a la Kermani and Plett 
            % 2001, eq. 9.
            %
            [D,L,R] = physics.getEigensystemAt(state1,state2);
            D12 = diag(D);
            D = abs(D12);
            % Fix:
            D1 = diag(physics.getEigensystemAt(state1));
            D2 = diag(physics.getEigensystemAt(state2));
            epsilon = max(0,max(D12-D1,D2-D12));
            i = D < epsilon;
            D(i) = epsilon(i);
            % AFC:
            D = R*(D.*L);
        end
        %% AFC diffusion: Roe + Hoffmann-Chiang (NOT RECOMMENDED)
        function D = diffusionRoeHoffmannChiang(state1,state2,physics)
            % Computes a tensorial numerical diffusion matrix block via AFC
            % a la Kuzmin et al. 2012, eq. 42; additionally, applies the
            % entropy fix of Hoffmaann and Chiang a la Kermani and Plett 
            % 2001, eq. 10.
            %
            [D,L,R] = physics.getEigensystemAt(state1,state2);
            D = abs(diag(D));
            % Fix:
            epsilon = .1; % usual range: 0 < epsilon < .125 (problem-dependent)
            i = D < epsilon;
            D(i) = (D(i).^2+epsilon.^2)./(2*epsilon);
            % AFC:
            D = R*(D.*L);
        end
        %% AFC diffusion: Roe + Kermani-Plett
        function D = diffusionRoeKermaniPlett(state1,state2,physics)
            % Computes a tensorial numerical diffusion matrix block via AFC
            % a la Kuzmin et al. 2012, eq. 42; additionally, applies the
            % entropy fix of Kermani and Plett 2001, eq. 11 (version 2).
            %
            [D,L,R] = physics.getEigensystemAt(state1,state2);
            D12 = diag(D);
            D = abs(D12);
            % Fix:
            D1 = diag(physics.getEigensystemAt(state1));
            D2 = diag(physics.getEigensystemAt(state2));
            epsilon = 4*max(0,max(D12-D1,D2-D12));
            i = D < epsilon;
            D(i) = (D(i).^2+epsilon(i).^2)./(2*epsilon(i));
            % AFC:
            D = R*(D.*L);
        end
        %% AFC diffusion a la Kuzmin et al. 2012, eq. 44
        function D = diffusionArithmetic(state1,state2,physics)
            [D,L,R] = physics.getEigensystemAt(.5*(state1+state2));
            D = R*abs(D)*L;
        end
        %% AFC diffusion a la Kuzmin et al. 2012, eq. 45
        function D = diffusionScalar(state1,state2,physics) % Rusanov-like
            D = physics.getEigensystemAt(state1,state2);
            D = max(abs(D(:)))*eye(size(D));
        end
        %% AFC diffusion a la Kuzmin et al. 2012, eq. 46
        function D = diffusionRobust(state1,state2,physics) % true Rusanov (i.e. local Lax-Friedrichs)
            d1 = physics.getEigensystemAt(state1);
            d1 = max(abs(d1(:)));
            d2 = physics.getEigensystemAt(state2);
            d2 = max(abs(d2(:)));
            D = max(d1,d2)*eye(physics.equationCount);
        end
    end
end