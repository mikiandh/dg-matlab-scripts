classdef WENO < Limiter_legendre
    % Simple WENO limiter for DG by Zhong and Zhu (2013), with additional 
    % improvements from Zhu et al. (2016).
    properties (Constant)
        % WENO linear weights determine the amount of up/down-winding done.
        % Their sum must add 1; the smaller the value for left and right
        % neighbors, the smaller the loss of accuracy in smooth regions,
        % but the weaker the limiting near discontinuities.
        linearWeights = [.001 .998 .001]
        % Small value to avoid division by zero.
        epsilon = 1e-6
    end
    properties (Access = protected)
        % Quadrature-based smoothness indicator requires an array of
        % DG (i.e. Legendre) basis instances (one per element), where 
        % useful quadrature and derivative data live.
        bases
        element2basis
    end
    methods
        %% Constructor
        function this = WENO(varargin)
            this@Limiter_legendre(varargin{:});
        end
        %% Apply (initialization, extension)
        function applyInitial(this,mesh,solver)
            % Instantiates a DG version of each basis in the mesh. Then
            % reverts to the superclass's default implementation.
            %
            % Preallocation:
            this.element2basis = ones(1,mesh.elementCount);
            this.bases = repelem(DG,1,numel(mesh.bases));
            % Fill arrays of unique bases and connectivity:
            for k = 1:numel(this.bases)
                this.bases(k) = DG(mesh.bases(k).basisCount);
                this.element2basis([mesh.elements.basis] == mesh.bases(k)) = k;
            end
            % Continue as default:
            this.applyInitial@Limiter_legendre(mesh,solver)
        end
    end
    methods (Access = protected)
        %% Apply (equation-wise)
        function applyNoSync(this)
            % Applies the WENO limiter asynchronously (i.e. on each
            % characteristic variable independently). Characteristic
            % decomposition has been carried out by the superclass.
            %
            % Compute modified Legendre coefficients of each troubled cell:
            this.coefsL(:,1,:) = this.coefs(:,1,:); % enforces a matching mean value, leaves the rest unchanged
            this.coefsR(:,1,:) = this.coefs(:,1,:);
            % Evaluate smoothnesses of each troubled cell:
            
        end
    end
    methods (Static, Access = protected)
        %% Smoothness indicators
        function [betaL,beta,betaR] = getSmoothness(I,p,basis,QL,Q,QR)
            % Evaluate smoothnesses of a troubled cell and its 2 neighbors:
            aux = 2.^(2*(1:p)'-1);
            aux = (QL*basis.derivatives)^2*basis.gaussWeights;
        end
    end
end