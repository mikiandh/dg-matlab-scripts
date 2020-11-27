function plotFun_single(inputCell)
figure(1)
[n,N,M,t,y0,y,z0,z,f,fy0,fy,fz0,fz] = inputCell{:};
subplot(1,2,1)
    plot(t,y0,t,y,t,z0,t,z)
    xlabel('Position')
    ylabel('Sample value')
    legend('q','q^h','\partial_t q','\partial_t q^h','Location','northoutside')
subplot(1,2,2)
    bar(f,[abs(fy0); abs(fy); abs(fz0); abs(fz)].'/M)
    hold on
    plot(n*[1 -1;1 -1],[0 0; 1 1],'--k')
    hold off
    legend('F(q)','F(q^h)','F(\partial_t q)','F(\partial_t q^h)','Target','Location','northoutside')
    xlabel('Wavemode')
    ylabel('Amplitude')
    xlim([-1 1]/2 + [0 N])
drawnow('limitrate')
end