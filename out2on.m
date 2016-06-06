function mon=out2on(mout,Lmax)
% mon=out2on(mout,Lmax)
%
% Transforms a list of coefficients in the formn of ADDMOUT into the form
% of ADDMON, that is from m=[0 -1 0 1 -2 -1 0 -1 2] into 
% m=[0 0 -1 1 0 -1 1 -2 2] such that it can be easily transformed into
% lmcosi.
%
% Last modified by plattner-at-princeton.edu, 07/25/2012


mon=nan(size(mout));

[~,~,~,~,mzero]=addmon(Lmax);
[~,~,mzout]=addmout(Lmax);


for l=0:Lmax
    mon(mzero(l+1)-1,:)=mout(mzout(l+1),:);
    for m=1:l
        mon(mzero(l+1)+(2*m-2),:)=mout(mzout(l+1)-m,:);
        mon(mzero(l+1)+(2*m-1),:)=mout(mzout(l+1)+m,:);
    end
end