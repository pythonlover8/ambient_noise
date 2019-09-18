% stack daily noise files
clear 
close all

sps=5;
load 'day_xcorrs/XCtest_3C_day_corrs_1.mat';

C_corrs_ZZ = zeros(size(day_corrs_ZZ));
C_corrs_EE = C_corrs_ZZ;
C_corrs_NN = C_corrs_ZZ;
C_corrs_EN = C_corrs_ZZ;
C_corrs_NE = C_corrs_ZZ;
N_SRpairs  = length(day_corrs_ZZ(:,1));
Ndays      = zeros(N_SRpairs,1);

%% 
nfiles=10;
for ifile = 1:nfiles
    load([ 'day_xcorrs/XCtest_3C_day_corrs_' num2str(ifile) '.mat']);
    
    for isrp = 1:N_SRpairs
        if sum(isnan(day_corrs_ZZ(isrp,:)))==0
            C_corrs_ZZ(isrp,:) = C_corrs_ZZ(isrp,:) + day_corrs_ZZ(isrp,:);
            Ndays(isrp) = Ndays(isrp)+1;
        end
    end
    
    for isrp = 1:N_SRpairs
        if sum(isnan(day_corrs_EE(isrp,:)))==0
            C_corrs_EE(isrp,:) = C_corrs_EE(isrp,:) + day_corrs_EE(isrp,:);
        end
    end
    
    for isrp = 1:N_SRpairs
        if sum(isnan(day_corrs_NN(isrp,:)))==0
            C_corrs_NN(isrp,:) = C_corrs_NN(isrp,:) + day_corrs_NN(isrp,:);
        end
    end
    
    for isrp = 1:N_SRpairs
        if sum(isnan(day_corrs_EN(isrp,:)))==0
            C_corrs_EN(isrp,:) = C_corrs_EN(isrp,:) + day_corrs_EN(isrp,:);
        end
    end
    
    for isrp = 1:N_SRpairs
        if sum(isnan(day_corrs_NE(isrp,:)))==0
            C_corrs_NE(isrp,:) = C_corrs_NE(isrp,:) + day_corrs_NE(isrp,:);
        end
    end
    
    file_frac = ifile/nfiles
end

for isrp = 1:N_SRpairs
    if Ndays(isrp)>3
        C_corrs_ZZ(isrp,:) = C_corrs_ZZ(isrp,:)./Ndays(isrp); 
        C_corrs_EE(isrp,:) = C_corrs_EE(isrp,:)./Ndays(isrp); 
        C_corrs_NN(isrp,:) = C_corrs_NN(isrp,:)./Ndays(isrp); 
        C_corrs_EN(isrp,:) = C_corrs_EN(isrp,:)./Ndays(isrp); 
        C_corrs_NE(isrp,:) = C_corrs_NE(isrp,:)./Ndays(isrp); 
    else
        C_corrs_ZZ(isrp,:) = NaN;
        C_corrs_EE(isrp,:) = NaN;
        C_corrs_NN(isrp,:) = NaN;
        C_corrs_EN(isrp,:) = NaN;
        C_corrs_NE(isrp,:) = NaN;
    end
end


%% filter only for plotting to check RR and TT results, usually comment out and save full bandwidth stacked correlations
%
sps=5;
dt=1/sps;
Wn=[0.04,0.16].*(dt*2); 
[b,a]=butter(1,Wn);  
r1=35; % secs
r=r1/(dt*.5*length(C_corrs_ZZ(1,:)));
%Tf=tukeywin(length(C_corrs_ZZ(1,:)),r);

for isrp = 1:length(C_corrs_ZZ(:,1)) 
    Tf=tukeywin(length(C_corrs_ZZ(isrp,~isnan(C_corrs_ZZ(isrp,:)))),r);
    C_corrs_ZZ(isrp,~isnan(C_corrs_ZZ(isrp,:))) = filtfilt(b,a,C_corrs_ZZ(isrp,~isnan(C_corrs_ZZ(isrp,:))).*Tf');
    %C_corrs_EE(isrp,:) = filtfilt(b,a,C_corrs_EE(isrp,:).*Tf');
    %C_corrs_NN(isrp,:) = filtfilt(b,a,C_corrs_NN(isrp,:).*Tf');
    %C_corrs_EN(isrp,:) = filtfilt(b,a,C_corrs_EN(isrp,:).*Tf');
    %C_corrs_NE(isrp,:) = filtfilt(b,a,C_corrs_NE(isrp,:).*Tf');
end
%}

%% plot ZZ, EE, NN, EN, NE
tt = -600:1/sps:600;

ex_index = 384; % just choosing one to plot

figure(1), clf
subplot(5,1,1)
plot(tt, C_corrs_ZZ(ex_index,:)./max(C_corrs_ZZ(ex_index,:)),'k-')
axis([-600 600 -1.1 1.1])
title('ZZ')
subplot(5,1,2)
plot(tt, C_corrs_EE(ex_index,:)./max(C_corrs_EE(ex_index,:)),'k-')
axis([-600 600 -1.1 1.1])
title('EE')
subplot(5,1,3)
plot(tt, C_corrs_NN(ex_index,:)./max(C_corrs_NN(ex_index,:)),'k-')
axis([-600 600 -1.1 1.1])
title('NN')
subplot(5,1,4)
plot(tt, C_corrs_EN(ex_index,:)./max(C_corrs_EN(ex_index,:)),'k-')
axis([-600 600 -1.1 1.1])
title('EN')
subplot(5,1,5)
plot(tt, C_corrs_NE(ex_index,:)./max(C_corrs_NE(ex_index,:)),'k-')
axis([-600 600 -1.1 1.1])
title('NE')

% rotate
for isrp=1:length(C_corrs_ZZ(:,1))
    [dist, theta] = distance(vslat(isrp),vslon(isrp),vrlat(isrp),vrlon(isrp));
    [dist, phi] = distance(vrlat(isrp),vrlon(isrp),vslat(isrp),vslon(isrp));
    
    TT(isrp,:) = -cosd(theta)*cosd(phi)*C_corrs_EE(isrp,:) + cosd(theta)*sind(phi)*C_corrs_EN(isrp,:) - sind(theta)*sind(phi)*C_corrs_NN(isrp,:) + sind(theta)*cosd(phi)*C_corrs_NE(isrp,:);
    RR(isrp,:) = -sind(theta)*sind(phi)*C_corrs_EE(isrp,:) - sind(theta)*cosd(phi)*C_corrs_EN(isrp,:) - cosd(theta)*cosd(phi)*C_corrs_NN(isrp,:) - cosd(theta)*sind(phi)*C_corrs_NE(isrp,:);
end

figure(2), clf
subplot(3,1,1)
plot(tt, C_corrs_ZZ(ex_index,:)./max(C_corrs_ZZ(ex_index,:)),'k-')
axis([-600 600 -1.1 1.1])
title('ZZ')
subplot(3,1,2)
plot(tt, RR(ex_index,:)./max(RR(ex_index,:)),'k-')
axis([-600 600 -1.1 1.1])
title('RR')
subplot(3,1,3)
plot(tt, TT(ex_index,:)./max(TT(ex_index,:)),'k-')
axis([-600 600 -1.1 1.1])
title('TT')
save(['stacked.mat'],C_corrs_ZZ)
%% plot all ZZ RR TT

[ds,di] = sort(sr_dist);

figure(3), clf
subplot(1,3,1)
imagesc(C_corrs_ZZ(di,:))
colormap('jet')
caxis([-0.05 0.05])

subplot(1,3,2)
imagesc(RR(di,:))
colormap('jet')
caxis([-0.05 0.05])

subplot(1,3,3)
imagesc(TT(di,:))
colormap('jet')
caxis([-0.05 0.05])


