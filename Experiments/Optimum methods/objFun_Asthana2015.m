function F = objFun_Asthana2015(basis)
% Objective function from Asthana & Jameson (2015).
%
% I/O
%  basis: discretization to evaluate
%  ---
%  F: cost <-> badness <-> penalty (smaller is better)
%
t = 100;

% %%% Vanilla version %%%
%     function f = fun(k)
%         [z,k,e] = basis.getFourierFootprint('criterion','energy','wavenumbers',k);
%         f = abs(1 - exp(1i*t*(k - 1i*z))).*e;
%     end
% %%%%%%%%%%%%%%%%%%%%%%%

%%% Combined mode version %%%
    function f = fun(k)
        [~,angs,amps] = basis.getAngAmp(t,'wavenumbers',k);
        f = abs(1 - amps.*exp(1i*angs));
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

F = integral(@fun,0,pi*basis.basisCount);
F = F/(basis.basisCount*pi); % I prefer this normalization
end