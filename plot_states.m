%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot state boundary
load us_states ;
n=length(state);
hold on ; 
for ct=1:n
  plot(state(ct).polygon(:,1),state(ct).polygon(:,2), 'color' ,[0.4 0.4 0.4] );
  hold on ;
end