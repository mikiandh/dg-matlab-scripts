classdef BDF < Limiter
    % Hierarchical moment limiter by Biswas, Devine and Flaherty (1994).
    % Free of user-defined parameters. Does not guarantee any TVD or TVB 
    % property (not even in the means).
    %
    properties (Access = protected, Constant)
        % Biswas et al., 1994 suggest carrying out 2 limiter sweeps each
        % time the limiter is called. This ensures that the modes used to
        % determine if limiting is required have themselves been limited.
        sweepCount = 2
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
            % Limiter sweeps:
            for sweep = 1:this.sweepCount
                % Get/set limited Legendre coefficients of each element:
                coefs{2} = this.getLimitedLegendre(elements(1));
                for k = 2:length(elements)-1
                    coefs{1} = this.getLimitedLegendre(elements(k));
                    elements(k-1).basis.setLegendre(elements(k-1),coefs{2});
                    coefs{2} = coefs{1};
                end
                elements(end).basis.setLegendre(elements(end),coefs{2});
            end
        end
    end
    methods (Access = protected)
        %% Limit an element
        function Q = getLimitedLegendre(this,element)
            % Applies hierarchical limiting a la Biswas et al., 1994, on an
            % element. All Legendre coefficients are limited in a single 
            % vectorized operation and then returned. The element's
            % solution is not modified, but the troubled components are
            % marked as limited.
            %
            % Aliases:
            I = this.physics.equationCount;
            N = element.dofCount;
            % Get unlimited legendre coefficents:
            Q = element.basis.getLegendre(element); % current element
            %%% TO DO %%%
        end
    end
    methods (Static)
        %% Information
        function info = getInfo
            info = 'BDF';
        end
    end
end