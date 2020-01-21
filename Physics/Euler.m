classdef Euler < Physics
    properties
        riemannSolver
    end
    methods
        %% Constructor
        function euler = Euler(boundaryConditions,riemannSolver)
            euler.equationCount = 3;
            euler.boundaryConditionsFunction = @Euler.applyTransmissiveBoundaryConditions;
            euler.riemannSolver = @Euler.riemannRoe;
            if nargin > 0
                switch boundaryConditions
                    case {'reflecting','reflective'}
                        euler.boundaryConditionsFunction = @euler.applyReflectingBoundaryConditions;
                    case 'transmissive'
                        euler.boundaryConditionsFunction = @euler.applyTransmissiveBoundaryConditions;
                    otherwise
                        error('Unknown boundary condition.')
                end
            end
            if nargin > 1
                switch riemannSolver
                    case 'exact'
                        euler.riemannSolver = @Euler.riemannExact;
                    case 'Rusanov'
                        euler.riemannSolver = @Euler.riemannRusanov;
                    case {'Roe','Roe1'}
                        euler.riemannSolver = @Euler.riemannRoe1;
                    case 'Roe0'
                        euler.riemannSolver = @Euler.riemannRoe0;
                    case 'HLL'
                        euler.riemannSolver = @Euler.riemannHLL;
                    case 'HLLE'
                        euler.riemannSolver = @Euler.riemannHLLE;
                    case 'HLLC'
                        euler.riemannSolver = @Euler.riemannHLLC;
                    otherwise
                        error('Unknown Riemann solver.')
                end
            end
        end
    end
    methods
        %% Apply boundary conditions
        function applyBoundaryConditions(this,mesh)
            this.boundaryConditionsFunction(mesh);
        end
        %% Numerical flux
        function [flux, waveSpeeds] = riemannFlux(this,stateL,stateR)
            [flux, waveSpeeds] = this.riemannSolver(stateL,stateR);
        end
        %% Reflecting boundary conditions
        function applyReflectingBoundaryConditions(this,mesh)
            % Right boundary:
            state = Euler.stateToPrimitive(mesh.elements(end).stateR);
            state(2) = -state(2);
            state = Euler.primitiveToState(state);
            [mesh.elements(end).riemannR,waveSpeeds] = this.riemannFlux(mesh.elements(end).stateR,state);
            mesh.edges{end}.computeTimeDeltas(waveSpeeds);
            % Left boundary:
            state = Euler.stateToPrimitive(mesh.elements(1).stateL);
            state(2) = -state(2);
            state = Euler.primitiveToState(state);
            [mesh.elements(1).riemannL,waveSpeeds] = this.riemannFlux(state,mesh.elements(1).stateL);
            mesh.elements(1).riemannL = -mesh.elements(1).riemannL;
            mesh.edges{1}.computeTimeDeltas(waveSpeeds);
        end
        %% Transmissive boundary conditions
        function applyTransmissiveBoundaryConditions(this,mesh)
            % Right-most edge:
            [mesh.elements(end).riemannR,waveSpeeds] = this.riemannFlux(...
                mesh.elements(end).stateR,mesh.elements(end).stateR);
            mesh.edges{end}.computeTimeDeltas(waveSpeeds);
            % Left-most edge:
            [mesh.elements(1).riemannL,waveSpeeds] = this.riemannFlux(...
                mesh.elements(1).stateL,mesh.elements(1).stateL);
            mesh.elements(1).riemannL = -mesh.elements(1).riemannL;
            mesh.edges{1}.computeTimeDeltas(waveSpeeds);
        end
    end
    methods (Static)
        %% Flux function
        function fluxes = flux(states)
            fluxes = states;
            states = Euler.stateToPrimitive(states);
            fluxes(1,:) = fluxes(2,:);
            fluxes(2,:) = fluxes(2,:).*states(2,:) + states(3,:);
            fluxes(3,:) = states(2,:).*(fluxes(3,:) + states(3,:));
        end
        %% Compute vector(s) of all primitive variables from state vector(s)
        function vars = allPrimitiveVarsFromState(state)
            vars = zeros(5,size(state,2));
            vars(1,:) = state(1,:); % density
            vars(2,:) = state(2,:)./vars(1,:); % velocity
            vars(6,:) = state(3,:)./vars(1,:); % total energy per unit mass
            vars(7,:) = 0.5*vars(2,:).^2; % kinetic energy per unit mass
            vars(8,:) = vars(6,:) - vars(7,:); % internal energy per unit mass
            vars(3,:) = 0.4*(state(3,:) - vars(1,:).*vars(7,:)); % pressure
            vars(4,:) = 1.4*vars(6,:) - 0.4*vars(7,:); % total enthalpy per unit mass (includes kinetic energy)
            vars(5,:) = sqrt(0.4*vars(4,:) - vars(7,:)); % speed of sound
        end
        %% Overwrite state vector(s) with vector(s) of (reduced) primitive variables
        function state = stateToPrimitive(state)
            state(2,:) = state(2,:)./state(1,:); % velocity
            state(3,:) = 0.4*(state(3,:) - 0.5*state(1,:).*state(2,:).^2); % pressure
        end
        %% Overwrite vector(s) of (reduced) primitive variables with state vector(s)
        function vars = primitiveToState(vars)
            vars(2,:) = vars(1,:).*vars(2,:); % momentum
            vars(3,:) = vars(3,:)/0.4 + 0.5*vars(2,:).^2./vars(1,:); % total energy
        end
        %% Primitive variables (one by one) from state vector(s)
        function [r,u,p,a,H] = getPrimitivesFromState(state)
            r = state(1,:);
            u = state(2,:)/state(1,:);
            p = 0.4*(state(3,:) - 0.5*r.*u.^2);
            a = realpow(1.4*p/r,0.5);
            H = (state(3,:) + p)/r;
        end
        %% Riemann solver (exact)
        function [flux,S] = riemannExact(stateL,stateR)
            % Wave speed estimate:
            stateL = Euler.stateToPrimitive(stateL);
            stateR = Euler.stateToPrimitive(stateR);
            S = max(abs(stateL(2)) + sqrt(1.4*stateL(3)/stateL(1)),...
                abs(stateR(2)) + sqrt(1.4*stateR(3)/stateR(1)));
            % Numerical flux:
            [r,u,p] = riemannEulerExact(1e-6,0,...
                stateL(1),stateL(2),stateL(3),...
                stateR(1),stateR(2),stateR(3),0);
            flux = Euler.flux(Euler.primitiveToState([r u p]'));
        end
        %% Riemann solver (Rusanov)
        function [flux,S] = riemannRusanov(stateL,stateR)
            this.flux1 = Euler.flux(stateL) + Euler.flux(stateR);
            this.flux2 = stateR - stateL;
            % Wave speed estimate:
            
            %%% TESTING %%%
            pL = 0.4*(stateL(3) - 0.5*stateL(2)^2/stateL(1));
            pR = 0.4*(stateR(3) - 0.5*stateR(2)^2/stateR(1));
            PASS = Euler.checkWithinPhysicalBounds(stateL,stateR);
            %%%%%%%%%%%%%%%
            
            [~,uL,~,aL,~] = Euler.getPrimitivesFromState(stateL); %%% FIXME: Woodward & Colella bursts here
            [~,uR,~,aR,~] = Euler.getPrimitivesFromState(stateR);
            S = max(abs(uL) + aL, abs(uR) + aR);
            % Numerical flux:
            flux = 0.5*(this.flux1 - S*this.flux2);
        end
        %% Riemann solver (Roe-Pike without entropy fix)
        function [flux,lambdas] = riemannRoe0(stateL,stateR)
            % Primitive variables and jumps:
            [rL,uL,pL,~,~] = Euler.getPrimitivesFromState(stateL);
            [rR,uR,pR,~,~] = Euler.getPrimitivesFromState(stateR);
            dr = rR - rL;
            du = uR - uL;
            dp = pR - pL;
            % Roe state vectors:
            roeL = stateL/realpow(rL,0.5);
            roeL(3) = 1.4*roeL(3) - 0.2*roeL(2)^2/roeL(1);
            roeR = stateR/realpow(rR,0.5);
            roeR(3) = 1.4*roeR(3) - 0.2*roeR(2)^2/roeR(1);
            % Roe averages:
            r = 1/(roeL(1) + roeR(1));
            u = r*(roeL(2) + roeR(2)); % velocity
            H = r*(roeL(3) + roeR(3)); % total enthalpy
            a = realpow(0.4*(H - 0.5*u^2),0.5); % speed of sound
            r = roeL(1)*roeR(1); % density
            % Eigenvalues:
            lambdas = [u-a; u; u+a];
            % Eigenvectors:
            K = [1 1 1; lambdas'; H-u*a 0.5*u^2 H+u*a];
            % Weights:
            alphas = [0.5*(dp-r*a*du)/a^2; dr - dp/a^2; 0.5*(dp+r*a*du)/a^2];
            % Numerical flux:
            flux = Euler.flux(stateL) + Euler.flux(stateR);
            flux = flux - K*(abs(lambdas).*alphas);
            flux = 0.5*flux;
        end
        %% Riemann solver (Roe-Pike with Harten-Hyman entropy fix)
        function [flux,lambdas] = riemannRoe1(stateL,stateR)
            % Primitive variables and jumps:
            [rL,uL,pL,aL,~] = Euler.getPrimitivesFromState(stateL);
            [rR,uR,pR,aR,~] = Euler.getPrimitivesFromState(stateR);
            dr = rR - rL;
            du = uR - uL;
            dp = pR - pL;
            % Roe state vectors:
            roeL = stateL/realpow(rL,0.5);
            roeL(3) = 1.4*roeL(3) - 0.2*roeL(2)^2/roeL(1);
            roeR = stateR/realpow(rR,0.5);
            roeR(3) = 1.4*roeR(3) - 0.2*roeR(2)^2/roeR(1);
            % Roe averages:
            r = 1/(roeL(1) + roeR(1));
            u = r*(roeL(2) + roeR(2)); % velocity
            H = r*(roeL(3) + roeR(3)); % total enthalpy
            a = realpow(0.4*(H - 0.5*u^2),0.5); % speed of sound
            r = roeL(1)*roeR(1); % density
            % Eigenvalues:
            lambdas = [u-a; u; u+a];
            % Eigenvectors:
            K = [1 1 1; lambdas'; H-u*a 0.5*u^2 H+u*a];
            % Weights:
            alphas = [0.5*(dp-r*a*du)/a^2; dr - dp/a^2; 0.5*(dp+r*a*du)/a^2];
%             del = stateR - stateL;
%             alphas = del;
%             alphas(2) = 0.4/a^2*(del(1)*(H - u^2) + u*del(2) - del(3));
%             alphas(1) = 0.5/a*(del(1)*(u + a) - del(2) - a*alphas(2));
%             alphas(3) = del(1) - alphas(1) - alphas(2);
            % Entropy fix (Harten-Hyman):
            aStar = rL + alphas(1);
            uStar = (stateL(2) + alphas(1)*(u - a))/aStar;
            pStar = 0.4*(stateL(3) + alphas(1)*(H - u*a) - 0.5*uStar^2*aStar);
            aStar = realpow(1.4*pStar/aStar,0.5);
            if uL < aL && uStar > aStar % left sonic rarefaction
                sL = uL - aL;
                sR = uStar - aStar;
                lambdas(1) = sL*(sR - lambdas(1))/(sR - sL);
            else % recalculate star region from the right
                aStar = rR - alphas(3);
                uStar = (stateR(2) - alphas(3)*(u + a))/aStar;
                pStar = 0.4*(stateR(3) - alphas(3)*(H + u*a) - 0.5*uStar^2*aStar);
                aStar = realpow(1.4*pStar/aStar,0.5);
                if -uR < aR && -uStar > aStar % right sonic rarefaction
                    sL = uStar + aStar;
                    sR = uR + aR;
                    lambdas(3) = sR*(lambdas(3) - sL)/(sR - sL);
                end
            end
            % Numerical flux:
            flux = Euler.flux(stateL) + Euler.flux(stateR);
            flux = flux - K*(abs(lambdas).*alphas);
            flux = 0.5*flux;
            %%%
            %i = lambdas < 0;
            %flux = Euler.flux(stateL) + K(:,i)*(lambdas(i).*alphas(i));
            %%%
            %i = lambdas > 0;
            %flux = Euler.flux(stateR) - K(:,i)*(lambdas(i).*alphas(i));
        end
        %% Riemann solver (HLL)
        function [flux,lambdas] = riemannHLL(stateL,stateR)
            % Primitive variables:
            [rL,uL,pL,aL,~] = Euler.getPrimitivesFromState(stateL);
            [rR,uR,pR,aR,~] = Euler.getPrimitivesFromState(stateR);
            % Wave speed estimates (PVRS-based, eqs. 10.59 - 10.62)
            p = 0.5*(pR + pL) + 0.125*(uL - uR)*(rL + rR)*(aL + aR);
            p = max(0,p);
            if p <= pL
                sL = uL - aL;
            else
                sL = realpow((1 + 6*p/pL)/7,0.5);
                sL = uL - aL*sL;
            end
            if p <= pR
                sR = uR + aR;
            else
                sR = realpow((1 + 6*p/pR)/7,0.5);
                sR = uR + aR*sR;
            end
            lambdas = [sL; sR];
            % Numerical flux (eq. 10.21):
            if sL >= 0
                flux = Euler.flux(stateL);
            elseif sR > 0
                flux = (sR*Euler.flux(stateL) - sL*Euler.flux(stateR) + ...
                    sL*sR*(stateR - stateL))/(sR - sL);
            else
                flux = Euler.flux(stateR);
            end
        end
         %% Riemann solver (HLLE, according to LeVeque)
        function [flux,lambdas] = riemannHLLE(stateL,stateR)
            % Primitive variables and jumps:
            [rL,uL,~,aL,~] = Euler.getPrimitivesFromState(stateL);
            [rR,uR,~,aR,~] = Euler.getPrimitivesFromState(stateR);
            % Roe state vectors:
            roeL = stateL/realpow(rL,0.5);
            roeL(3) = 1.4*roeL(3) - 0.2*roeL(2)^2/roeL(1);
            roeR = stateR/realpow(rR,0.5);
            roeR(3) = 1.4*roeR(3) - 0.2*roeR(2)^2/roeR(1);
            % Roe averages:
            r = 1/(roeL(1) + roeR(1));
            u = r*(roeL(2) + roeR(2)); % velocity
            H = r*(roeL(3) + roeR(3)); % total enthalpy
            a = realpow(0.4*(H - 0.5*u^2),0.5); % speed of sound
            % Wave speed estimates (LeVeque, eq. 15.60)
            sL = min(uL - aL,u - a);
            sR = max(uR + aR,u + a);
            lambdas = [sL; sR];
            % Numerical flux (eq. 10.21):
            if sL >= 0
                flux = Euler.flux(stateL);
            elseif sR > 0
                flux = (sR*Euler.flux(stateL) - sL*Euler.flux(stateR) + ...
                    sL*sR*(stateR - stateL))/(sR - sL);
            else
                flux = Euler.flux(stateR);
            end
        end
        %% Riemann solver (HLLC)
        function [flux,lambdas] = riemannHLLC(stateL,stateR)
            % Primitive variables:
            [rL,uL,pL,aL,~] = Euler.getPrimitivesFromState(stateL); %%% FIXME: WoodwardColella bursts at t = 0.027 here
            [rR,uR,pR,aR,~] = Euler.getPrimitivesFromState(stateR);
            % Wave speed estimates (eqs. 10.67 - 10.70)
            p = 0.5*(pR + pL) + 0.125*(uL - uR)*(rL + rR)*(aL + aR);
            p = max(0,p);
            if p <= pL
                sL = uL - aL;
            else
                sL = realpow((1 + 6*p/pL)/7,0.5);
                sL = uL - aL*sL;
            end
            if p <= pR
                sR = uR + aR;
            else
                sR = realpow((1 + 6*p/pR)/7,0.5);
                sR = uR + aR*sR;
            end
            sC =  pR - pL + rL*uL*(sL - uL) - rR*uR*(sR - uR);
            sC = sC/(rL*(sL - uL) - rR*(sR - uR));
            lambdas = [sL; sC; sR];
            % Numerical flux (eqs. 10.71 and 10.75 - 10.76)
            if sL >= 0
                flux = Euler.flux(stateL);
            elseif sC > 0
                flux = pL + rL*(sL - uL)*(sC - uL) + ...
                    pR + rR*(sR - uR)*(sC - uR);
                flux = 0.5*sL*flux*[0; 1; sC];
                flux = sC*(sL*stateL - Euler.flux(stateL)) + flux;
                flux = flux/(sL - sC);
            elseif sR > 0
                flux = pL + rL*(sL - uL)*(sC - uL) + ...
                    pR + rR*(sR - uR)*(sC - uR);
                flux = 0.5*sR*flux*[0; 1; sC];
                flux = sC*(sR*stateR - Euler.flux(stateR)) + flux;
                flux = flux/(sR - sC);
            else
                flux = Euler.flux(stateR);
            end
        end
        %% Return 3D arrays of left and right eigenvectors
        function [L,R] = getEigenvectors(meanStates)
            % Returns the eigenvector matrices corresponding to each
            % mean state in the input, including left-averaged and
            % right-averaged (left) eigenvectors (Roe average).
            %
            % meanStates(i,j,k)
            % L.{'left','none','right'}(i,j,k)
            % R.{'left','none','right'}(i,j,k)
            %  L struct fields: left-averaged, not averaged, right-averaged
            %  i: system dimension
            %  j: eigenvector index
            %  k: element
            % Preallocate:
            R.left = ones(3,3,size(meanStates,3)-2); % 1st and last states are excluded
            %%%R.left = repmat(eye(3,3),1,1,size(meanStates,3)-2);
            R.none = R.left;
            R.right = R.left;
            L.left = R.left;
            L.none = R.left;
            L.right = R.left;
            %%%return %%% TESTING
            for k = 2:size(meanStates,3)-1
                % Not-averaged eigenvectors:
                [r,u,~,a,H] = Euler.getPrimitivesFromState(meanStates(:,1,k));
                R.none(2:3,:,k-1) = [u-a u u+a; H-u*a 0.5*u^2 H+u*a];
                L.none(:,:,k-1) = 0.2/a^2*[H+2.5*a*(u-a) -(u+2.5*a) 1; -2*H+10*a^2 2*u -2; H-2.5*a*(u+a) -u+2.5*a 1];
                % Left-averaged eigenvectors (Roe average):
                r = realpow(r,0.5);
                [rAvg,uAvg,~,~,HAvg] = Euler.getPrimitivesFromState(meanStates(:,1,k-1));
                rAvg = realpow(rAvg,0.5);
                aux = 1/(rAvg + r);
                uAvg = (rAvg*uAvg + r*u)*aux;
                HAvg = (rAvg*HAvg + r*H)*aux;
                aAvg = realpow(0.4*(HAvg - 0.5*uAvg^2),0.5);
                R.left(2:3,:,k-1) = [uAvg-aAvg uAvg uAvg+aAvg; HAvg-uAvg*aAvg 0.5*uAvg^2 HAvg+uAvg*aAvg];
                L.left(:,:,k-1) = 0.2/aAvg^2*[HAvg+2.5*aAvg*(uAvg-aAvg) -uAvg-2.5*aAvg 1; -2*HAvg+10*aAvg^2 2*uAvg -2; HAvg-2.5*aAvg*(uAvg+aAvg) -uAvg+2.5*aAvg 1];
                % Right-averaged eigenvectors (Roe average):
                [rAvg,uAvg,~,~,HAvg] = Euler.getPrimitivesFromState(meanStates(:,1,k+1));
                rAvg = realpow(rAvg,0.5);
                aux = 1/(rAvg + r);
                uAvg = (rAvg*uAvg + r*u)*aux;
                HAvg = (rAvg*HAvg + r*H)*aux;
                aAvg = realpow(0.4*(HAvg - 0.5*uAvg^2),0.5);
                R.right(2:3,:,k-1) = [uAvg-aAvg uAvg uAvg+aAvg; HAvg-uAvg*aAvg 0.5*uAvg^2 HAvg+uAvg*aAvg];
                L.right(:,:,k-1) = 0.2/aAvg^2*[HAvg+2.5*aAvg*(uAvg-aAvg) -uAvg-2.5*aAvg 1; -2*HAvg+10*aAvg^2 2*uAvg -2; HAvg-2.5*aAvg*(uAvg+aAvg) -uAvg+2.5*aAvg 1];
            end
        end
        %% Return left and right eigenvectormatrices for given mean states
        function [L,R] = getEigenvectorsAt(meanStates)
            % Returns the eigenvector matrices evaluated at the mean state
            % passed as input.
            %
            % meanStates(i,j,k)
            % L(i,j,k)
            % R(i,j,k)
            %  i: system dimension
            %  j: eigenvector index
            %  k: element
            %
            K = size(meanStates,3);
            R = zeros(3,3,K);
            L = R;
            for k = 1:K
                [~,u,~,a,H] = Euler.getPrimitivesFromState(meanStates(:,1,k));
                R(:,:,k) = [1 1 1; u-a u u+a; H-u*a 0.5*u^2 H+u*a];
                L(:,:,k) = 0.2/a^2*[H+2.5*a*(u-a) -(u+2.5*a) 1; -2*H+10*a^2 2*u -2; H-2.5*a*(u+a) -u+2.5*a 1];
            end
        end
        %% Return left and right eigenvectormatrices for given mean states
        function [L,R] = getLocalEigenvectors(stateL,stateR)
            % Returns the eigenvector matrices evaluated at the local
            % Roe-averaged state between the passed mean states.
            %
            % stateL(i,1), stateR(i,1)
            % L(i,j)
            % R(i,j)
            %  i: system dimension
            %  j: eigenvector index
            %
            % Primitive variables:
            [r,u,~,~,H] = Euler.getPrimitivesFromState(stateL);
            [rAvg,uAvg,~,~,HAvg] = Euler.getPrimitivesFromState(stateR);
            % Roe averages:
            r = realpow(r,0.5);
            rAvg = realpow(rAvg,0.5);
            aux = 1/(rAvg + r);
            uAvg = (rAvg*uAvg + r*u)*aux;
            HAvg = (rAvg*HAvg + r*H)*aux;
            aAvg = realpow(0.4*(HAvg - 0.5*uAvg^2),0.5);
            % Local (a.k.a averaged) eigenvectors:
            R = [uAvg-aAvg uAvg uAvg+aAvg; HAvg-uAvg*aAvg 0.5*uAvg^2 HAvg+uAvg*aAvg];
            L = 0.2/aAvg^2*[HAvg+2.5*aAvg*(uAvg-aAvg) -uAvg-2.5*aAvg 1; -2*HAvg+10*aAvg^2 2*uAvg -2; HAvg-2.5*aAvg*(uAvg+aAvg) -uAvg+2.5*aAvg 1];
        end
        %% Check for positivity in certain quantities
        function PASS = checkWithinPhysicalBounds(stateL,stateR)
            PASS = 0;
            try
                [rL,uL,pL,aL,~] = Euler.getPrimitivesFromState(stateL);
                [rR,uR,pR,aR,~] = Euler.getPrimitivesFromState(stateR);
            catch
                return
            end
            if rL <= 0 || rR <= 0 % negative density
                return
            elseif pL <= 0 || pR <= 0 % negative pressure
                return
            elseif stateL(3) <= 0 || stateR(3) <= 0 % negative total energy
                return
            elseif 5*(aL + aR) <= uR - uL % vacuum generation (Toro, eq. 4.40)
                return
            end
            PASS = 1;
        end
    end
end