function varargout=plotplates(c11,cmn,res)
% [handl,XY]=PLOTPLATES(c11,cmn,res)
% 
% Plots plate boundaries contained within 'c11' and 'cmn'
%
% res   : 0 regular (default) 
%       : 1 regular (for compatibility with PLOTCONT)
%       : 2 Global Mollweide projection
%
% Last modified by fjsimons-at-alum.mit.edu, 24.11.2004

defval('c11',[0 90])
defval('cmn',[360 -90])
defval('res',0)

switch res
 case {0,1}
  pathname=getenv('IFILES');
  fid=fopen(fullfile(pathname,'PLATES','plates.mtl'),'r','b');
  plates=fread(fid,[1692 2 ],'uint16');
  plates(plates==0)=NaN;
  plates=plates/100-90;
  lon=plates(:,1);
  lat=plates(:,2);
  lon(~(lon>=c11(1) & lon<=cmn(1)))=NaN;
  lat(~(lat<=c11(2) & lat>=cmn(2)))=NaN;
 case 2
  load(fullfile(getenv('IFILES'),'COASTS','platm'))
  lon=plxm;
  lat=plym;
end

hold on

handl=plot(lon,lat,'k');
XY=[lon(:) lat(:)];

varnames={'handl' 'XY'};
for index=1:nargout
  varargout{index}=eval(varnames{index});
end

%file=fullfile(getenv('IFILES'),'PLATES','nuvel1_plates.datp');
%[err,dat]=unix([ 'awk ''{if($1~/>/) {print "NaN NaN"} else {print}}'' ',...
%      file]);
%dat=str2num(dat);
%[m,n]=size(dat);
% Look at  [(1:m)' dat] find out where NaN's need be inserted.
%positions=[24 28 44 106 410];
%lat=insert(dat(:,1),NaN,positions);
%lon=insert(dat(:,2),NaN,positions);
%plot(lon,lat)
%plates=([lon' lat']+90)*100;
%pathname=getenv('IFILES');
%fid=fopen([pathname '/PLATES/plates.mtl'],'wb');
%fwrite(fid,plates,'uint16')
%fclose(fid)

