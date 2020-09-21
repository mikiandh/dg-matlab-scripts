function plotFun_single(inputCell)
[j,k,w,f] = inputCell{:};
yyaxis left
plot(k,real(w),':','DisplayName',sprintf('f = %g',f),'Tag',['disp_' num2str(j)])
yyaxis right
plot(k,imag(w),':','Tag',['diss_' num2str(j)])
drawnow('limitrate')
end