classdef BDF < Limiter_legendre
    % Hierarchical moment limiter by Biswas, Devine and Flaherty (1994).
    % Free of user-defined parameters. Does not guarantee any TVD or TVB 
    % property (not even in the means).
    methods
        %% Constructor
        function this = BDF(varargin)
            % Superclass constructor:
            this = this@Limiter_legendre(varargin{:});
        end
    end
    methods (Access = protected)
        %% Apply (equation-wise)
        function applyNoSync(this)
            % Applies the minmod limiter asynchronously (i.e. on each
            % characteristic variable separately).
            %
            % Loop over characteristic variables:
            for i = 1:this.I
                % Flag every element for limiting:
                k = true(this.K,1);
                % Limit hierarchically (top to bottom):
                for j = this.J:-1:1
                    % Get limited coefs.:
                    coefs_min = this.coefs(k,j+1,i);
                    this.coefs(k,j+1,i) = this.minmod(...
                        (2*j-1).*this.coefs(k,j+1,i),...
                        this.coefs(k,j,i) - this.coefsL(k,j,i),...
                        this.coefsR(k,j,i) - this.coefs(k,j,i))./(2*j-1);
                    % Update list of elements flagged for limiting:
                    k(k) = abs(this.coefs(k,j+1,i) - coefs_min) > 1e-10 | ~coefs_min;
                end
            end
        end
    end
    methods (Static)
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