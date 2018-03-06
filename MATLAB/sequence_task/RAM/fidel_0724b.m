% incidence of trisynapic response to EC stim as a function of stim
% frequency in Fidel
i2 = [0 0 0 0 0];
t2 = 0:1000/2:1000/2*(length(i2)-1);
i4 = [0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1];
t4 = 0:1000/4:1000/4*(length(i4)-1);
i6 = [0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1];
t6 = 0:1000/6:1000/6*(length(i6)-1);
i8 = [0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1];
t8 = 0:1000/8:1000/8*(length(i8)-1);
i10 = [0 0 0 0 0 0 0 0 0 1 0 1 1 1 1 1 1 1 1 0 1];
t10 = 0:1000/10:1000/10*(length(i10)-1);
i12 = [0 0 0 0 0 0 0 0 0 0 0 0 1 0 1 1 0 1 1 0 1 0 1 0 1 0 0 1 0 1];
t12 = 0:1000/12:1000/12*(length(i12)-1);
i14 = [0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 0 1 0 0 1 0 0 1 0 0 1 0 1];
t14 = 0:1000/14:1000/14*(length(i14)-1);
i16 = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 1 0 0 1 0 0 1 0 1 0 0 1 0];
t16 = 0:1000/16:1000/16*(length(i16)-1);

figure('Name','trisynaptic','NumberTitle','off','Units','normalized','Position',[1/4 1/4 1/2 1/2],'Color','w');
ind = find(i2); plot(t2(ind),ones(1,length(ind)),'ko','MarkerSize',12,'MarkerFaceColor','k'); hold on;
ind = find(~i2); plot(t2(ind),ones(1,length(ind)),'ko','MarkerSize',12);
ind = find(i4); plot(t4(ind),2*ones(1,length(ind)),'ko','MarkerSize',12,'MarkerFaceColor','k'); hold on;
ind = find(~i4); plot(t4(ind),2*ones(1,length(ind)),'ko','MarkerSize',12);
ind = find(i6); plot(t6(ind),3*ones(1,length(ind)),'ko','MarkerSize',12,'MarkerFaceColor','k'); hold on;
ind = find(~i6); plot(t6(ind),3*ones(1,length(ind)),'ko','MarkerSize',12);
ind = find(i8); plot(t8(ind),4*ones(1,length(ind)),'ko','MarkerSize',12,'MarkerFaceColor','k'); hold on;
ind = find(~i8); plot(t8(ind),4*ones(1,length(ind)),'ko','MarkerSize',12);
ind = find(i10); plot(t10(ind),5*ones(1,length(ind)),'ko','MarkerSize',12,'MarkerFaceColor','k'); hold on;
ind = find(~i10); plot(t10(ind),5*ones(1,length(ind)),'ko','MarkerSize',12);
ind = find(i12); plot(t12(ind),6*ones(1,length(ind)),'ko','MarkerSize',12,'MarkerFaceColor','k'); hold on;
ind = find(~i12); plot(t12(ind),6*ones(1,length(ind)),'ko','MarkerSize',12);
ind = find(i14); plot(t14(ind),7*ones(1,length(ind)),'ko','MarkerSize',12,'MarkerFaceColor','k'); hold on;
ind = find(~i14); plot(t14(ind),7*ones(1,length(ind)),'ko','MarkerSize',12);
ind = find(i16); plot(t16(ind),8*ones(1,length(ind)),'ko','MarkerSize',12,'MarkerFaceColor','k'); hold on;
ind = find(~i16); plot(t16(ind),8*ones(1,length(ind)),'ko','MarkerSize',12);
set(gca,'Box','off','XLim',[-100 2000],'YLim',[0 9],'YTick',1:8,'YTickLabel','2|4|6|8|10|12|14|16');
xlabel('time (ms)'); ylabel('stimulus frequency (Hz)');