function stop = plotFun_dispDiss(eta,optimValues,state,p,i,export)
% Plots the dispersion and dissipation relations(blue lines), as well as
% their ratio (a la Adams et al., 2015), in red, and the cutoff wavenumber
% (a la Moura et al. 2015), in black, of the current iteration in an FR
% optimization problem. Compares them with those of DG (dotted lines).
stop = false;
switch state
    case 'init' % set up the plot
        
        [r,k,kM] = FR('DG',p).getDispDissRatios('combined');
        kc = k(find(imag(kM) <= -.01,1,'first')); % cheap estimation (good enough)
        
        hold on
        yyaxis left
        ylabel '\Re(\kappa_m), \Im(\kappa_m)'
        plot([kc kc],[min(imag(kM)) max(real(kM))],':k')
        plot([kc kc],[min(imag(kM)) max(real(kM))],'-k')
        
        plot([0 k(end)],[0 k(end)],'--')
        plot(k,real(kM),':')
        plot(k,real(kM),'-')
        
        plot([0 k(end)],[0 0],'--')
        plot(k,imag(kM),':')
        plot(k,imag(kM),'-')
        
        yyaxis right
        ylabel '|d\Re(\kappa_m)/d\kappa - 1|/|\Re(\kappa_m)|'
        plot(k,r,':')
        plot(k,r,'-')
        
        hold off
        xlabel '\kappa'
    case 'iter'
        basis = FR({'eta',eta},p);
        [r,k,kM] = basis.getDispDissRatios('combined');
        kc = k(find(imag(kM) <= -.01,1,'first'));
        
        title(sprintf('%s - iter.: %d',basis.getName,optimValues.iteration))
        yyaxis left
        h = get(gca,'Children');
        yyaxis right
        h = [h; get(gca,'Children')];
        set(h([1 4 7 9]),{'XData','YData'},{
            k,imag(kM)
            k,real(kM)
            [kc kc],[min(imag(kM)) max(real(kM))]
            k,r
            })
        if export.gif
            try
                im = frame2im(getframe(gcf)); % will fail inside a parfor
            catch
                return
            end
            [imind,cm] = rgb2ind(im,256);
            if optimValues.iteration == 0
                imwrite(imind,cm,[export.name '_' num2str(i) '.gif'],'gif','Loopcount',inf);
            else
                imwrite(imind,cm,[export.name '_' num2str(i) '.gif'],'gif','WriteMode','append');
            end
        end
end
end