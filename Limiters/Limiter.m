classdef Limiter < handle
    properties (Abstract)
        everyStage
    end
    methods (Abstract)
        apply(this,mesh,timeDelta)
    end
    methods (Static)
        %% Factories
        % No limiter (placeholder)
        function limiter = none(varargin)
            limiter = [];
        end
        % TVBM limiter, Cockburn & Shu 1989
        function limiter = TVBM(physics,M)
            if isa(physics,'Advection') || isa(physics,'Burgers')
                limiter = TVBM_scalar(M);
            elseif isa(physics,'Wave') || isa(physics,'Euler')
                limiter = TVBM_vector(physics,M);
            else
                error('Limiter not available.')
            end
        end
        % Moment limiter, Biswas et al 1994
        function limiter = Biswas(physics)
            if isa(physics,'Advection') || isa(physics,'Burgers')
                limiter = Biswas_scalar;
            elseif isa(physics,'Wave') || isa(physics,'Euler')
                limiter = Biswas_vector(physics);
            else
                error('Limiter not available.')
            end
        end
        % Moment limiter, Burbeau et al 2001
        function limiter = Burbeau(physics)
            if isa(physics,'Advection') || isa(physics,'Burgers')
                limiter = Burbeau_scalar;
            elseif isa(physics,'Wave') || isa(physics,'Euler')
                limiter = Burbeau_vector(physics);
            else
                error('Limiter not available.')
            end
        end
        % Moment limiter, Krivodonova 2007
        function limiter = Krivodonova(physics,ratio)
            if nargin == 1
                ratio = 0;
            end
            if isa(physics,'Advection') || isa(physics,'Burgers')
                limiter = Krivodonova_scalar(ratio);
            elseif isa(physics,'Wave') || isa(physics,'Euler')
                limiter = Krivodonova_vector(physics,ratio);
            else
                error('Limiter not available.')
            end
        end
        % AP-TVD detector + PFGM limiter, Wang 2009
        function limiter = Wang(physics)
            if isa(physics,'Advection') || isa(physics,'Burgers')
                limiter = Wang_scalar;
            elseif isa(physics,'Wave') || isa(physics,'Euler')
                limiter = Wang_vector(physics);
            else
                error('Limiter not available.')
            end
        end
        % Algebraic Flux Correction, Kuzmin et. al. 2012
        function limiter = AFC(physics)
            if isa(physics,'Advection') || isa(physics,'Burgers')
                limiter = AFC_scalar(physics);
            elseif isa(physics,'Wave') || isa(physics,'Euler')
                warning('Vector version of AFC is still under development.')
                limiter = AFC_scalar(physics); % TO DO: implement AFC_vector
            else
                error('Limiter not available.')
            end
        end
        %% Minmod function for 3 scalar arguments
        function d = minmod(a,b,c)
            d = sign(a);
            if d == sign(b) && d == sign(c)
                d = d*min(abs([a b c]));
            else
                d = 0;
            end
        end
        %% Minmod function for 3 arguments which are 1D column vectors
        function d = minmod1D(a,b,c)
            d = zeros(length(a),1);
            M = abs(horzcat(a,b,c));
            a = sign(a);
            ids = a == sign(b) & a == sign(c);
            d(ids) = min(M(ids,:),[],2);
            d = a.*d;
        end
        %% Modified minmod function of Cockburn & Shu (1D arrays)
        function d = modminmod1D(m,a,b,c)
            d = zeros(length(a),1);
            M = abs(horzcat(a,b,c));
            a = sign(a);
            ids1 = M(:,1) <= m;
            d(ids1) = M(ids1,1);
            ids2 = a == sign(b) & a == sign(c) & ~ids1;
            d(ids2) = min(M(ids2,:),[],2);
            d = a.*d;
        end
        %% Minmod function for 3 arguments that are 2D arrays
        function D = minmod2D(A,B,C)
            [I,J] = size(A);
            A = reshape(A,1,[]);
            B = reshape(B,1,[]);
            C = reshape(C,1,[]);
            M = vertcat(A,B,C);
            A = sign(A);
            B = sign(B);
            C = sign(C);
            M = abs(M);
            ids = A == B & A == C;
            D = zeros(1,J*I);
            D(:,ids) = min(M(:,ids));
            D = A.*D;
            D = reshape(D,I,J);
        end
    end
end