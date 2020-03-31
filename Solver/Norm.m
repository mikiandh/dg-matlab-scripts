classdef Norm < handle
    % This class defines various finite-dimensional vector space norms, 
    % pseudo-norms, and the like. Norms are normalized by the domain length
    % (a la Wang et al. 2013, e.g. equation 10).
    %
    properties (Constant)
        tol = 1e-3 % tolerance for iterative norms
    end
    properties (SetAccess = protected)
        p
        subs
        vals
        rows
    end
    enumeration
        L1 (1)
        L2 (2)
        Max (inf)
        ErrorL1 (1)
        ErrorL2 (2)
        ErrorMax (inf)
        Mass (1)
        TVM (nan)
        TV (nan)
    end
    methods
        %% Constructor
        function this = Norm(p)
            this.p = p;
        end
        %% Evaluate
        function compute(this,mesh,fun)
            % Computes norms of each scalar component of the solution on
            % a given mesh (according to the norm type). Supports being 
            % called on an array of norm instances. Returns 'true' iff the
            % norm computation was successful.
            %
            % Arguments
            %  mesh: grid where the approximate solution "lives"
            %  fun: exact solution (function of space only)
            %
            % Update norm properties:
            [this.rows] = deal(size(mesh.elements(1).states,1));
            [this.vals] = deal(nan(this(1).rows,1));
            [this.subs] = deal(reshape(mesh.getBreakLocations(false),2,[])); % polynomial endpoints
            % Loop over norm instances:
            for norm = this
                % De facto polymorphism:
                switch norm
                    case {Norm.L1, Norm.L2}
                        % Solution p-norms:
                        norm.pNorm(@(x) mesh.sample(x))
                    case {Norm.Max}
                        % Solution maximum norm (estimate):
                        norm.maxNorm(@(x) mesh.sample(x))
                    case {Norm.ErrorL1, Norm.ErrorL2}
                        % Error p-norms:
                        norm.pNorm(@(x) mesh.sample(x) - fun(x))
                    case {Norm.ErrorMax}
                        % Error maximum norm:
                        norm.maxNorm(@(x) mesh.sample(x) - fun(x))
                    case {Norm.Mass}
                        % Solution mass pseudo-norm:
                        norm.massNorm(mesh)
                    case {Norm.TV}
                        % Solution total variation (estimate):
                        norm.tvNorm(mesh)
                    case {Norm.TVM}
                        % Solution total variation in the means:
                        norm.tvmNorm(mesh)
                end
            end
        end
    end
    methods (Access = protected)
        %% Mass norms
        function massNorm(this,mesh)
            % Computes the "mass norm" (non-normalized L1, without absolute
            % values) of solution vector components via Gauss quadrature of
            % optimal order at each break span.
            %
            % Computation of patch averages exploits the fact that, for
            % bases that satify the partition of unity property at nodes
            % (e.g. Legendre, BSpline), each lumped mass matrix component
            % is equal to the integral of its associated basis function 
            % over the whole patch.
            %
            % For Legendre, Lagrange and BSpline bases, results are exact.
            %
            % Preallocate:
            this.vals = zeros(this.rows,1);
            % Loop over elements:
            for element = mesh.elements
                % "Polymorphism":
                if isa(element.basis,'Legendre') % cell-averages directly available
                    this.vals = this.vals + element.states(:,1)*element.dx;
                elseif isa(element.basis,'Lagrange') || isa(element.basis,'BSpline') % nodal partition of unity
                    this.vals = this.vals + element.states*sum(element.basis.massMatrix,2)*.5*element.dx;
                else % general case (approximate)
                    this.vals = this.vals + element.getLegendre(1)*element.dx;
                end
            end
        end
        %% L^p norms
        function pNorm(this,fun)
            % Computes a p-norm of the given function using this norm's
            % properties. Applied subinterval-wise.
            this.vals = Algorithms.quadvgk(@(x) abs(fun(x)).^this.p,this.subs,this.rows)/(this.subs(end) - this.subs(1));
            this.vals = nthroot(this.vals,this.p);
        end
        %% Maximum norm
        function maxNorm(this,fun)
            % Estimates the maximum norm of a given function by 
            % finding its global maximum via global search and local
            % minimization. Applied mesh-wide.
            %
            % Instantiate a global search:
            gs = GlobalSearch('Display','off','MaxTime',1,'FunctionTolerance',this.tol);
            % Instantiate a local minimization problem:
            problem = createOptimProblem('fmincon',...
                'x0',median(this.subs(:)),'lb',this.subs(1),'ub',this.subs(end),...
                'options',optimoptions('fmincon','Algorithm','active-set','FunctionTolerance',this.tol));
            % Loop over function components:
            for i = 1:this.rows
                % Set the (scalar) function to minimize:
                fun_i = @(y) y(i);
                problem.objective = @(x) -abs(fun_i(fun(x)));
                % Solve the global minimization:
                [~,this.vals(i)] = run(gs,problem);
            end
            % Retrieve maxima:
            this.vals = -this.vals;
        end
        %% Total variation in the means
        function tvmNorm(this,mesh)
            % Computes, equation-wise, the total variation in the means of
            % an approximate solution (i.e. total variation a la FV, see 
            % Leveque, 2002, p. 109, applied to 1st Legendre coefficients).
            %
            % Loop over edges (so that ghost elements are also included):
            this.vals(1:this.rows,1) = 0;
            for edge = mesh.edges
                % Approximate Legendre coefficients patch-wide:
                QL = edge.elementL.getLegendre(1,1:this.rows);
                QR = edge.elementR.getLegendre(1,1:this.rows);
                % Accumulate TV in Legendre coefficients:
                this.vals = this.vals + abs(QR - QL);
            end
        end
        %% TV estimate
        function tvNorm(this,mesh)
            % Compute the total variation of the current approximate
            % solution stored in the mesh. Convergence acceleration using
            % linear Richardson extrapolation.
            %
            maxIters = 30;
            x = mesh.getGaussLocations; % global quadrature point locations
            q = mesh.sample(x); % state vector at every quadrature location (initial samples)
            function tv = TV
                tv = sum(abs(diff(q,1,2)),2);
            end
            A(:,2) = TV;
            R0 = inf;
            % Iterate until reaching tolerance or excessive cost:
            for iter = 1:maxIters
                [q,x] = resampleBisection(q,x);
                A = horzcat(A(:,2),TV);
                if all(abs(diff(A,1,2)) < 1e-10)
                    this.vals = A(:,end);
                    return
                end
                this.vals = Algorithms.richardsonExtrapolate(A);
                if all(abs(this.vals - R0) < 1e-10)
                    return
                end
                if length(x) > 1e6 || any(isnan(this.vals)) || any(this.vals > 1e200)
                    %warning('Failed to converge on tolerance.')
                    break
                end
                R0 = this.vals;
            end
            function [q,x] = resampleBisection(q0,x0)
                % Refine (in a nested way) a set of solution samples (from
                % this mesh), by adding to the provided sampled states and
                % locations new ones obtained via bisection. Only the new
                % samples are computed.
                %
                x = zeros(1,2*length(x0)-1);
                q = zeros(size(q0,1),length(x));
                ids = logical(mod(1:length(x),2)); % old entries
                x(ids) = x0;
                q(:,ids) = q0;
                ids = ~ids; % new entries
                x(ids) = .5*(x0(2:end) + x0(1:end-1));
                q(:,ids) = mesh.sample(x(ids));
            end
        end
    end
end