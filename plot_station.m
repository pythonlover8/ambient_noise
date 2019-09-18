clear all
iday=2;
file_name = ['/Users/a123/Desktop/PHD/amb_noi/test_3C_xcorr/day_xcorrs/XCtest_3C_day_corrs_' num2str(iday) '.mat'];
load(file_name)
file_name1=['/Users/a123/Desktop/PHD/amb_noi/3C_xcorr/ZZ_day_corrs_' num2str(iday) '.mat'];
load(file_name1)
load us_states;
n=length(state);
index=[11,17,31,52,61,121,151,406];
close all;
for ii = 1:length(index)
    h(ii)=figure;
    subplot(2,2,1:2);
    for ct=1:n
    plot(state(ct).polygon(:,1),state(ct).polygon(:,2), 'color' ,[0.4 0.4 0.4] );
    hold on ;
    end
    plot(vslon(index(ii)),vslat(index(ii)),'r.','MarkerSize',20);
    plot(vrlon(index(ii)),vrlat(index(ii)),'r.','MarkerSize',20);
    text(vslon(index(ii))+1,vslat(index(ii))+1,stn_list(vsrc_num(index(ii))));
    text(vrlon(index(ii))+1,vrlat(index(ii))+1,stn_list(vrec_num(index(ii))));
    xlim([-120,-100])
    ylim([30,50]);
    subplot(2,2,3)
    plot(data(index(ii),:));
    subplot(2,2,4)
    plot(day_corrs_ZZ(index(ii),:));
    filename="/Users/a123/Desktop/PHD/amb_noi/3C_xcorr/record_day"
    filename=strcat(filename,num2str(iday),'_',num2str(index(ii)),".png");
    saveas(h(ii),filename);
end