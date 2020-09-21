function stop = plotFun_eigenvalues(p,eta,optimValues,state,varargin)
% Plots the fourier footprint of the current iterant of an FR optimization
% problem.
stop = false;
switch state
    case 'init' % set up the plot
        z = FR('DG',p).getFourierFootprint;
        set(gca,{'XAxisLocation','YAxisLocation'},{'origin','origin'})
        plot(z(:),'.k','DisplayName','eta = 0');
        hold on
        title(sprintf('Degree: %d',p))
        xlabel \Re(z*)
        ylabel \Im(z*)
        legend('-DynamicLegend','Location','SouthWest')
    case 'iter'
        z = FR({'eta',eta},p).getFourierFootprint;
        if optimValues.iteration > 0
            delete(findobj(get(gca,'Children'),'Tag','old'));
            hNew = findobj(get(gca,'Children'),'Tag','new');
            hOld = copy(hNew,gca);
            set(hOld,{'Color','Tag'},{'b','old'});
            delete(hNew)
        end
        plot(z(:),'.r','DisplayName',sprintf('eta = %g (f = %g)',eta,optimValues.fval),'Tag','new');
end
end