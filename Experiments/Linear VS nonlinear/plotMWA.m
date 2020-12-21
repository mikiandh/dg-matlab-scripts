function plotMWA(inputCell)
figure(1)
[n,N,M,K,t,y0,y,z0,z,f,fy0,fy,fz0,fz] = inputCell{:};
subplot(2,2,1)
    plot(t,y0,t,y)
    ylabel('Solution')
    xlabel('Position')
    legend('q','q^h')
subplot(2,2,2)
    plot(t,z0,t,z)
    ylabel('Residual')
    xlabel('Position')
    legend('\partial_t q','\partial_t q^h')
subplot(2,2,3)
    bar(f,[abs(fy0); abs(fy)].'/M)
    hold on
    plot(n*[1 -1;1 -1],[0 0; 1 1],'--k')
    hold off
    legend('F(q)','F(q^h)','Target')
    xlabel('Wavemode')
    ylabel('Solution amplitude')
    xlim([-1 1]/2 + [0 N])
    nmults = unique([-n - 0:K:N -n + 0:K:N n - 0:K:N n + 0:K:N]);
    nmults(nmults < 0) = [];
    title(['n = [' num2str(nmults) ']'])
subplot(2,2,4)
    bar(f,[abs(fz0); abs(fz)].'/M)
    hold on
    plot(n*[1 -1;1 -1],[0 0; 1 1],'--k')
    hold off
    legend('F(\partial_t q)','F(\partial_t q^h)','Target')
    xlabel('Wavemode')
    ylabel('Residual amplitude')
    xlim([-1 1]/2 + [0 N])
drawnow('limitrate')
end