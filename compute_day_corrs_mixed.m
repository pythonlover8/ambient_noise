function [ day_corrs ] = compute_day_corrs_mixed( day_seis1, day_seis2, samprate, max_lag, vsrc_num, vrec_num )

% compute correlations
seis_npts = length(day_seis1(1,:));
seg_min_dur = 239;
seg_npts   = seg_min_dur*60*samprate;
seg_starts = 1:round(0.5*seg_min_dur*60*samprate):seis_npts-seg_npts-1;
nsegs      = length(seg_starts);
day_corrs  = zeros(length(vsrc_num),2*max_lag+1,nsegs);

for iseg = 1:nsegs
    seg_seis1 = day_seis1(:,seg_starts(iseg):seg_starts(iseg)+seg_npts-1); 
    seg_seis2 = day_seis2(:,seg_starts(iseg):seg_starts(iseg)+seg_npts-1); 
   
    for istn = 1:length(day_seis1(:,1)) % whiten non-NaN segs
        if sum(isnan(seg_seis1(istn,:)))==0 && sum(isnan(seg_seis2(istn,:)))==0
            seg_seis1(istn,:) = seg_seis1(istn,:) - mean(seg_seis1(istn,:));
            seg_seis1(istn,:) = detrend(seg_seis1(istn,:));
            seg_seis1(istn,:) = whiten(seg_seis1(istn,:),0.98,samprate); 
            seg_seis1(istn,:) = seg_seis1(istn,:).*tukeywin(seg_npts,0.01)';
            
            seg_seis2(istn,:) = seg_seis2(istn,:) - mean(seg_seis2(istn,:));
            seg_seis2(istn,:) = detrend(seg_seis2(istn,:));
            seg_seis2(istn,:) = whiten(seg_seis2(istn,:),0.98,samprate); 
            seg_seis2(istn,:) = seg_seis2(istn,:).*tukeywin(seg_npts,0.01)';
            
            %{
            clf
            plot(seg_seis(istn,:))
            %}
        end
    end
 
    seg_corrs = zeros(length(vsrc_num),2*max_lag+1);
    for isrp = 1:length(vsrc_num) %par 
        temp = xcorr(seg_seis1(vsrc_num(isrp),:),seg_seis2(vrec_num(isrp),:),max_lag); 
        temp = temp./max(abs(temp));
        seg_corrs(isrp,:) = temp;
    end
    
    day_corrs(:,:,iseg) = seg_corrs;
end

day_corrs = nanmean(day_corrs,3);

end

