function [ day_corrs ] = compute_day_corrs( day_seis, samprate, max_lag, vsrc_num, vrec_num )

% compute correlations
seis_npts = length(day_seis(1,:));
seg_min_dur = 239;
seg_npts   = seg_min_dur*60*samprate;
seg_starts = 1:round(0.5*seg_min_dur*60*samprate):seis_npts-seg_npts-1;
nsegs      = length(seg_starts);
day_corrs  = zeros(length(vsrc_num),2*max_lag+1,nsegs);

for iseg = 1:nsegs
    seg_seis = day_seis(:,seg_starts(iseg):seg_starts(iseg)+seg_npts-1); 
    
    for istn = 1:length(day_seis(:,1)) % whiten non-NaN segs
        if sum(isnan(seg_seis(istn,:)))==0
            seg_seis(istn,:) = seg_seis(istn,:) - mean(seg_seis(istn,:));
            seg_seis(istn,:) = detrend(seg_seis(istn,:));
            seg_seis(istn,:) = whiten(seg_seis(istn,:),0.98,samprate); 
            seg_seis(istn,:) = seg_seis(istn,:).*tukeywin(seg_npts,0.01)';
            
            %{
            clf
            plot(seg_seis(istn,:))
            %}
        end
    end
 
    seg_corrs = zeros(length(vsrc_num),2*max_lag+1);
    for isrp = 1:length(vsrc_num) %par 
        temp = xcorr(seg_seis(vsrc_num(isrp),:),seg_seis(vrec_num(isrp),:),max_lag); 
        temp = temp./max(abs(temp));
        seg_corrs(isrp,:) = temp;
    end
    
    day_corrs(:,:,iseg) = seg_corrs;
    
end

day_corrs = nanmean(day_corrs,3);

end

