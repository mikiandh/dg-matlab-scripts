classdef WENO < Limiter_legendre
    % Simple WENO limiter for DG by Zhong and Zhu (2013), with additional 
    % improvements from Zhu et al. (2016).
    properties (Constant)
        % WENO linear weights determine the amount of up/down-winding done.
        % Their sum must add 1; the smaller the value for left and right
        % neighbors, the smaller the loss of accuracy in smooth regions,
        % but the weaker the limiting near discontinuities.
        linearWeights = [.001 .998 .001]
    end
    methods
        %% Constructor
        function this = WENO(varargin)
            this@Limiter_legendre(varargin{:});
        end
    end
    methods (Access = protected)
        %% Apply (equation-wise)
        function applyNoSync(this)
            % Applies the WENO limiter asynchronously (i.e. on each
            % characteristic variable independently). Characteristic
            % decomposition is handled by the superclass.
            %
            % Compute modified Legendre coefficients of each troubled cell:
            this.coefsL(:,1,:) = this.coefs(:,1,:); % enforces a matching mean value, leaves the rest unchanged
            this.coefsR(:,1,:) = this.coefs(:,1,:);
            % Evaluate smoothnesses of every troubled cell and its 2 neighbors:
            w = this.getSmoothnesses(this.coefsL,this.coefs,this.coefsR);
            % Evaluate non-linear weights:
            w = this.linearWeights./(1e-6 + w).^2;
            w = w./sum(w,2);
            % Reconstruct WENO-limited Legendre coefficients:
            this.coefs = this.coefsL.*w(:,1,:) + this.coefs.*w(:,2,:) + this.coefsR.*w(:,3,:);
        end
    end
    methods (Static,Access = protected)
        %% Compute smoothness
        function betas = getSmoothnesses(varargin)
            % Evaluate smoothness indicator values for each given array of
            % Legendre coefficients.
            %
            % Arguments
            %  varargin: list of Legendre coefficient arrays (row: cell;
            %            column: basis function index; page: equation)
            % Output
            %  betas: smoothness indicator values (row: troubled cell;
            %         column: input index; page: equation)
            %
            % Deduce sizes (#troubled elements, #DoFs, #equations):
            [K,J,I] = size(varargin{1});
            % Preallocate smoothness indicator values:
            betas = zeros(K,nargin,I);
            % Precompute some stuff:
            [xi,w] = Legendre.quadratureGaussLegendre(J-1);
            phi = Legendre.getDerivatives(J,xi',J-1);
            % Loop over system components:
            for i = 1:I
                % Loop over derivatives (skipping the zeroth one):
                for s = 2:size(phi,3)
                    % Loop over cells within stencil (left, middle, right):
                    for l = 1:size(betas,2)
                        betas(:,l,i) = betas(:,l,i) + 2^(2*s-3).*(varargin{l}(:,:,i)*phi(:,:,s)).^2*w;
                    end
                end
            end
        end
    end
end