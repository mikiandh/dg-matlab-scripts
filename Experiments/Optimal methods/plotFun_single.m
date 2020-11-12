function plotFun_single(inputCell)
[j,k,w,name] = inputCell{:};
yyaxis left
plot(k,real(w),':','DisplayName',[name ' (disp.)'],'Tag',['disp_' num2str(j)])
yyaxis right
plot(k,imag(w),':','DisplayName',[name ' (diss.)'],'Tag',['diss_' num2str(j)])
drawnow('limitrate')
end