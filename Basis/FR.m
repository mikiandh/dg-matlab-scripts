classdef FR < Lagrange
    properties
        param % value for the "c" parameter in VCJH-class correction functions
        correctionsL % array of p+1 "gradient of left correction function" values
        correctionsR
        auxMassMatrix % mass matrix if FR was using Lagrange test functions (it is not)
    end
    methods
        %% Constructor
        function this = FR(param,varargin)
            this@Lagrange(varargin{:});
            this.auxMassMatrix = this.massMatrix;
            switch nargin
                case 0
                    this.param = 0; % defaults to DG if no correction function is specified
                    return;
                case 1
                    this.param = param;
                    return;
                otherwise
                    this.param = param;
                    if this.order == 1 % trivial case (1st order FV)
                        this.correctionsR = 0.5;
                        this.correctionsL = -0.5;
                        return
                    end
                % Derivative of the correction function at solution points:
                eta = this.getEtaParameter(param,this.order);
                dPn = Legendre.getLegendreAndDerivatives(this.degree+2,this.nodeCoords');
                this.correctionsR = this.rightVCJH(eta,dPn(this.degree:end,:));
                this.correctionsL = - flip(this.correctionsR);
                % Mass and gradient matrices (Dirac delta test functions):
                this.massMatrix = eye(this.basisCount);
                this.gradientMatrix = this.derivatives;
            end
        end
        %% Instantiate from prototype
        function fr = clone(this,p)
            fr = FR(this.param,p);
        end
        %% Oneliner info (extension)
        function name = getName(this)
            % Adds the VCJH correction function parameter to this basis's
            % one line descriptor.
            name = getName@Lagrange(this);
            aux = sprintf('%s(%s)',class(this),num2str(this.param));
            name = strrep(name,class(this),aux);
        end
        %% L2 projection (extension)
        function project(this,element,fun0)
            % Initializes an L2-preserving interpolant for FR, taking into
            % account the fact that its test functions are collocated Dirac
            % Deltas.
            %
            % Parent implementation:
            project@Lagrange(this,element,fun0);
            % Recover L2-preserving nodal values:
            element.states = element.states/this.auxMassMatrix;
        end
    end
    methods (Access = {?Basis,?Element})
        %% FR operator
        function computeResiduals(this,element,physics)
            element.computeFluxesFromStates(physics);
            element.interpolateFluxAtEdges;
            element.residuals = element.fluxes*this.derivatives +...
                (-element.riemannL - element.fluxL)*this.correctionsL +... % negative left Riemann flux to avoid the conventional sign assigned in the edge-based calculation
                (element.riemannR - element.fluxR)*this.correctionsR;
            element.residuals = -2/element.dx*element.residuals;
        end
    end
    methods (Static)
        %% Right VCJH function
        function g = rightVCJH(eta,legendre)
            % Evaluates the expression for the right correction function
            % given in Vincent et. al. 2011 (eq. 3.47) using provided parameter 
            % 'eta' and samples of the Legendre polynomials.
            %
            % If Legendre derivatives are given instead, output is the
            % derivative of the correction function, directly.
            %
            % Arguments
            %  eta: a constant parameter
            %  legendre: 2D array with 3 rows, matching Legendre
            %            polynomials p to p+2, where p+1 is the
            %            degree of the correction function; colums: sample
            %            locations.
            g = 0.5*(legendre(2,:) + (eta*legendre(1,:) + legendre(3,:))/(eta+1));
        end
        %% VCJH parameter 'eta' from degree and parameter 'c'
        function eta = getEtaParameter(param,degree)
            % Evaluates eq. 3.45 from Vincent et. al. 2011.
            %
            % Arguments
            %  param: VCJH parameter (can also be a character array)
            %  degree: degree of the correction function (NOT discretization)
            %
            p = factorial(degree);
            a = factorial(2*degree)/(2^degree*(p)^2);
            % VCJH-type scheme:
            if isa(param,'double')
                eta = 0.5*param*(2*degree+1)*(a*p)^2;
                return;
            end
            % Huyn-type schemes:
            switch param
                case 'min'
                    eta = -.5; % 'safe' minimum (halfway towards unsafe)
                case 'DG'
                    eta = 0;
                case 'Ga'
                    eta = degree/(degree+1);
                case 'LumpLo'
                    eta = (degree+1)/degree;
                case 'max'
                    eta = 1e16; % not-to-large number
                otherwise
                    error('Correction function unknown.')
            end
        end
        %% VCJH parameter 'c' from degree and parameter 'eta'
        function param = getParameter(param,degree)
            % Returns the value for 'c' that would result in a given 'eta'.
            %
            % Arguments
            %  param: VCJH parameter (in: eta; out: c)
            %  degree: degree of the correction function (NOT discretization)
            %
            p = factorial(degree);
            a = factorial(2*degree)/(2^degree*(p)^2);
            param = param/(0.5*(2*degree+1)*(a*p)^2);
        end
        %% Correction functions and derivatives
        function [gL,dgL,gR,dgR] = getCorrectionFunctionAndDerivative(c,p,xi)
            % Returns the p+1 degree correction functions and their 
            % derivative (left and right), obtained using a given "c" value
            % (VCJH-type correction functions), sampled at given locations. 
            %
            % Arguments
            %  c: VCJH parameter (e.g. c = 0 <-> DG)
            %  p: degree of the solution polynomial
            %  xi: sample locations (1D row array)
            %
            % Return
            %  gL,gR: 1D row arrays of correction function samples
            %  dgL,dgR: idem for derivatives
            %
            xi = xi(:)'; % ensure sample locations are a row vector
            if p == 0 % trivial case
                gL = -xi;
                dgL = -0.5+0*xi;
                gR = xi;
                dgR = 0.5+0*xi;
                return;
            end
            eta = FR.getEtaParameter(c,p);
            [dPn,Pn] = Legendre.getLegendreAndDerivatives(p+2,xi);
            Pn = Pn(p:end,:);
            gR = FR.rightVCJH(eta,Pn);
            gL = flip(gR);
            dPn = dPn(p:end,:);
            dgR = FR.rightVCJH(eta,dPn);
            dgL = - flip(dgR);
        end
    end
    methods (Access = protected)
        %% Assemble Fourier residual matrices (override)
        function [A,B,C] = getFourierMatrices(this,beta)
            % FR version.
            A = -2*this.gradientMatrix.' + sparse((1+beta)*this.correctionsL.'*this.left.' + (1-beta)*this.correctionsR.'*this.right.');
            B = -(1+beta)*sparse(this.correctionsL.'*this.right.'); % upwind
            C = -(1-beta)*sparse(this.correctionsR.'*this.left.'); % downwind
        end
    end
end