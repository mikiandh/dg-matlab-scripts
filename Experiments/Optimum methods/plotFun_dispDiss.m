function stop = plotFun_dispDiss(p,eta,optimValues,state,varargin)
% Plots the dispersion and dissipation relations of the current iteration
% in an FR optimization problem.
% Compares them with those of the previous iteration, DG's, and exact ones.
stop = false;
switch state
    case 'init' % set up the plot
        [z,k] = FR('DG',p).getFourierFootprint;
        hold on
        yyaxis left
        plot(k,k,'--k','DisplayName','Exact','Tag','exactDisp')
        plot(k,-imag(z(1,:)),':k','DisplayName','eta = 0','Tag','refDisp')
        ylabel \Re(\omega*)
        yyaxis right
        plot(k,0*k,'--k','Tag','exactDiss')
        plot(k,real(z(1,:)),':k','Tag','refDiss')
        ylabel \Im(\omega*)
        xlabel \kappa*
        xlim([0 inf])
        title(sprintf('Degree: %d',p))
    case 'iter'
        [z,k] = FR({'eta',eta},p).getFourierFootprint;
        if optimValues.iteration > 0
            for side = ["left" "right"]
                yyaxis(side)
                delete(findobj(get(gca,'Children'),'Tag','old'));
                hNew = findobj(get(gca,'Children'),'Tag','new');
                hOld = copy(hNew,gca);
                set(hOld,{'Color','Tag'},{'b','old'});
                delete(hNew)
            end
        end
        yyaxis left
        plot(k,-imag(z(1,:)),'-r',...
            'DisplayName',sprintf('eta = %g (f = %g)',eta,optimValues.fval),...
            'Tag','new')
        yyaxis right
        plot(k,real(z(1,:)),'-r','Tag','new')
        yyaxis left
        legend(get(gca,'Children'),'Location','South')
end
end