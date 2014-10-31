function sdwallconmovie
% SWALLCONMOVIE

L=18;
for index=1:361
  fnpl=sprintf('%s/SDWALLCONS-%i-%i.mat',...
	       fullfile(getenv('IFILES'),'WIECZOREK'),L,index);
  load(fnpl)                                             
  F(F<max(F(:))/100)=NaN;
  imagefnan([0 90],[360 -90],F,'default',minmax(F))
  % surf(F) ; shading flat; view(-180,70) ; axis off
  set(gca,'ytick',[-90:30:90])
  set(gca,'xtick',[0:90:360])
  deggies(gca)
  plotcont
  longticks(gca,2)
  a=title(sprintf('N = %i %s %i',0,'\rightarrow',index));
  set(a,'FontS',20)
  print('-depsc',...
	sprintf('/home/fjsimons/EPS/WREGCUM/ContCum_%3.3i',index))
end
