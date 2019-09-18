%% vertical noise xcorr
clear
close all

javaaddpath('IRIS-WS-2.0.6.jar');
addpath('day_data')
addpath('day_xcorrs')
sps=5;

% first run, get desired station list, then save for possible future runs
% Start = '2000-01-01';
% End =  '2001-12-31'; %one day before half were removed 
% stations = irisFetch.Stations('','XC','','','BHZ','StartTime',Start,'EndTime',End, ...
%                          'MinimumLatitude', 40 , 'MaximumLatitude', 50, 'MinimumLongitude', -113,'MaximumLongitude', -107);
% save('stations_XC','stations')
load stations_XC.mat
               
figure(1), clf
scatter([stations.Longitude],[stations.Latitude],'k^','filled')  
axis([ min([stations.Longitude]) max([stations.Longitude]) min([stations.Latitude]) max([stations.Latitude]) ])
hold on
plot_states

%%
for i=1:length(stations)
    stn_list{i} = stations(i).StationCode;
    net_list{i} = stations(i).NetworkCode;
    stn_lat(i)  = [stations(i).Latitude];
    stn_lon(i)  = [stations(i).Longitude];
    stn_elv(i)  = [stations(i).Elevation];
end
                   
n=1;
for ii = 1:length(stations)   
    for jj = ii:length(stations)
        if jj~= ii 
        vsrc_num(n) = ii;
        vslat(n)    = stn_lat(ii);
        vslon(n)    = stn_lon(ii);
        vselv(n)    = stn_elv(ii);
        
        vrec_num(n) = jj;
        vrlat(n)    = stn_lat(jj);
        vrlon(n)    = stn_lon(jj);
        vrelv(n)    = stn_elv(jj);
        n = n+1;
        end
    end
end

% sr dist 
sr_dist = zeros(size(vsrc_num));
for isrp = 1:length(vsrc_num)
    sr_dist(isrp) = distance(vslat(isrp),vslon(isrp),vrlat(isrp),vrlon(isrp));
end

% start times 
ndays = 10; % *****only using 10 days for test case
start_times = zeros(ndays,1);
start_times(1) = datestr2epoch('2000-10-01 00:00:00.0');

for ii = 2:length(start_times)
   start_times(ii) = start_times(ii-1) + 3600*24;
end

dt=1/sps;
day_t_vec  = 0:dt:24*3600; 

%%
max_lag = 600*sps;
req_dur  = 24*3600; % durection in seconds, to check for short data segments
for iday = 1:length(start_times)
    tic
    sprintf([ '\n Starting day ' num2str(iday)]) 
    
    % req or load
    %
    day_data_Z = NaN*zeros(length(stn_list),length(day_t_vec));
    day_data_E = day_data_Z;
    day_data_N = day_data_Z;
    
    for istn = 1:length(stn_lat)
        stn_name  = stn_list(istn); 
        netw_name = net_list(istn);
        st_time = epoch2datestr(start_times(iday));
        end_time = epoch2datestr(start_times(iday)+24*3600);
        traceZ = irisFetch.Traces(netw_name,stn_name,'*','BHZ,HHZ',st_time,end_time);
        traceE = irisFetch.Traces(netw_name,stn_name,'*','BHE,HHE',st_time,end_time);
        traceN = irisFetch.Traces(netw_name,stn_name,'*','BHN,HHN',st_time,end_time);
         
        % Z
        if length(traceZ)>=1
            if traceZ(1).sampleCount >= 0.5*req_dur*traceZ(1).sampleRate   
                isamprate = round(traceZ(1).sampleRate);
                dataZ = traceZ(1).data;     
                dataZ = dataZ - mean(dataZ);
                dataZ = detrend(dataZ); 
                idt=1/isamprate;
                Wn=[0.02,1].*(idt*2); 
                [b,a]=butter(1,Wn);  
                r1=50; 
                r=r1/(idt*.5*length(dataZ));
                Tf=tukeywin(length(dataZ),r);
                dataZ = filtfilt(b,a,dataZ.*Tf);  
                dataZ = resample(dataZ,sps,isamprate);   
                % interp onto standard time vector
                tr_st_time = traceZ(1).startTime;
                tr_st_time = datevec(tr_st_time);
                tr_st_sec  = 3600*tr_st_time(4) + 60*tr_st_time(5) + tr_st_time(6);
                tr_t_vec   = tr_st_sec:1/sps:24*3600+10;
                tr_t_vec   = tr_t_vec(1:length(dataZ));
                temp_data = interp1(tr_t_vec,dataZ,day_t_vec);
                day_data_Z(istn,:) = temp_data;
            end
        end
        
        % E
        if length(traceE)>=1
            if traceE(1).sampleCount >= 0.5*req_dur*traceE(1).sampleRate   
                isamprate = round(traceE(1).sampleRate);
                dataE = traceE(1).data;     
                dataE = dataE - mean(dataE);
                dataE = detrend(dataE); 
                idt=1/isamprate;
                Wn=[0.02,1].*(idt*2); %?%
                [b,a]=butter(1,Wn);  %?%
                r1=50; %?%
                r=r1/(idt*.5*length(dataE));
                Tf=tukeywin(length(dataE),r);
                dataE = filtfilt(b,a,dataE.*Tf);  
                dataE = resample(dataE,sps,isamprate);   
                % interp onto standard time vector
                tr_st_time = traceE(1).startTime;
                tr_st_time = datevec(tr_st_time);
                tr_st_sec  = 3600*tr_st_time(4) + 60*tr_st_time(5) + tr_st_time(6);
                tr_t_vec   = tr_st_sec:1/sps:24*3600+10;
                tr_t_vec   = tr_t_vec(1:length(dataE));
                temp_data = interp1(tr_t_vec,dataE,day_t_vec);
                day_data_E(istn,:) = temp_data;
            end
        end
        
        % N
        if length(traceN)>=1
            if traceN(1).sampleCount >= 0.5*req_dur*traceN(1).sampleRate    
                isamprate = round(traceN(1).sampleRate);
                dataN = traceN(1).data;     
                dataN = dataN - mean(dataN);
                dataN = detrend(dataN); 
                idt=1/isamprate;
                Wn=[0.02,1].*(idt*2); 
                [b,a]=butter(1,Wn);  
                r1=50; 
                r=r1/(idt*.5*length(dataN));
                Tf=tukeywin(length(dataN),r);
                dataN = filtfilt(b,a,dataN.*Tf);  
                dataN = resample(dataN,sps,isamprate);   
                % interp onto standard time vector
                tr_st_time = traceN(1).startTime;
                tr_st_time = datevec(tr_st_time);
                tr_st_sec  = 3600*tr_st_time(4) + 60*tr_st_time(5) + tr_st_time(6);
                tr_t_vec   = tr_st_sec:1/sps:24*3600+10;
                tr_t_vec   = tr_t_vec(1:length(dataN));
                temp_data = interp1(tr_t_vec,dataN,day_t_vec);
                day_data_N(istn,:) = temp_data;
            end
        end

        %{
        clf
        plot(day_t_vec,day_data_Z(istn,:),'k-')
        hold on
        plot(day_t_vec,day_data_E(istn,:),'g--')
        plot(day_t_vec,day_data_N(istn,:),'r--')
        keyboard
        %}

        sprintf([ '\n ' cell2mat(stn_list(istn)) ])
    end
    NaN_frac_Z = sum(isnan(day_data_Z(:)))/length(day_data_Z(:)) 
    NaN_frac_E = sum(isnan(day_data_E(:)))/length(day_data_E(:)) 
    NaN_frac_N = sum(isnan(day_data_N(:)))/length(day_data_N(:)) 
    
    save([ 'day_data/day_3C_data' num2str(iday) '.mat'],'day_data_Z','day_data_E','day_data_N','stn_list','net_list','stn_lat','stn_lon','stn_elv');
    %}
    
    load([ 'day_data/day_3C_data' num2str(iday) '.mat']);
    
    [ day_corrs_ZZ ] = compute_day_corrs( day_data_Z, sps, max_lag, vsrc_num, vrec_num );
    [ day_corrs_EE ] = compute_day_corrs( day_data_E, sps, max_lag, vsrc_num, vrec_num );
    [ day_corrs_NN ] = compute_day_corrs( day_data_N, sps, max_lag, vsrc_num, vrec_num );
    
    [ day_corrs_EN ] = compute_day_corrs_mixed( day_data_E, day_data_N, sps, max_lag, vsrc_num, vrec_num );
    [ day_corrs_NE ] = compute_day_corrs_mixed( day_data_N, day_data_E, sps, max_lag, vsrc_num, vrec_num );

    % save daily stacks
    file_name = ['day_xcorrs/XCtest_3C_day_corrs_' num2str(iday) '.mat'];
    save(file_name,'day_corrs_ZZ','day_corrs_EE','day_corrs_NN','day_corrs_EN','day_corrs_NE',...
                    'sr_dist','vslat','vslon','vselv','vsrc_num','vrlat','vrlon','vrelv','vrec_num','stn_list','net_list','-v7.3')
        
    toc
end
    
    
