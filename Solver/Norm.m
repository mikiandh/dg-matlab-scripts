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
                        norm.vals = mesh.getSolutionMass(1:norm.rows);
                    case {Norm.TV}
                        % Solution total variation (estimate):
                        norm.vals = mesh.getTotalVariation;
                    case {Norm.TVM}
                        % Solution total variation in the means:
                        norm.vals = mesh.getTVM(1:norm.rows);
                    otherwise
                        
                end
            end
        end
    end
    methods (Access = protected)
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
    end
end