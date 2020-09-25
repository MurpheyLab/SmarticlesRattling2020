% Simpler version of the symmetrize_smarticles function
% just finding a vector of rotation- and permutation- invariant observables

if(stRes==1); stIx=1; else; stIx=2; end %if re-analyzing, then store from start

%% Symmetrize data (For rotationally symmetric)
%relative coordinates to account for rotation and translation invariance:
relDat=repmat(mean(crdDat(Nwall+1:end,1:3,:)),[Nsm-Nwall,1,1])-crdDat(Nwall+1:end,1:3,:); %c.o.m. in each sm-cle's frame
% relDat=[crdDat(2:end,1:3,:);crdDat(1,1:3,:)]-crdDat(:,1:3,:);
% relDat1=[crdDat(end,1:3,:);crdDat(1:end-1,1:3,:)]-crdDat(:,1:3,:);
relDat(:,3,:)=mod(relDat(:,3,:)+pi,2*pi)-pi;
for si=1:Nsm-Nwall %rotate crd into the reference frames of prev sm.
    cth=cos(crdDat(Nwall+si,3,:)); sth=sin(crdDat(Nwall+si,3,:));
    for ti=1:size(relDat,3)
        relDat(si,1:2,ti)=relDat(si,1:2,ti)*[[cth(:,:,ti), -sth(:,:,ti)]; [sth(:,:,ti),cth(:,:,ti)]];
%         relDat1(si,1:2,ti)=relDat1(si,1:2,ti)*[[cth(:,:,ti), -sth(:,:,ti)]; [sth(:,:,ti),cth(:,:,ti)]];
    end
    if(size(crdDat,3)>1E6); si
    end
end %problem: translational d.o.f. are dominated by angular

%% Permutation-invarnat coordinates piDat
if(false) %store piDat
clear piDat;
% relDat=relDat1; relDat(:,2:3,:)=exp(1i.*relDat(:,2:3,:));
% relDat=[relDat; relDat1];
piDat(1,:,:)=mean(relDat(:,1:2,:));%-mean(mean(relDat),3);
for im=2:Nsm %find the first central moments
  piDat(im,:,:)=(moment(relDat(:,1:2,:),im)); piDat(im,:,:)=sign(piDat(im,:,:)).*abs(piDat(im,:,:)).^(1/im);
  if(size(crdDat,3)>1E6); im
  end
end
piDat=reshape(piDat,2*Nsm,[]); % coordinates in 6D space ~(<x>,<x^2>,<x^3>,<y>,<y^2>,<y^3>)
else %just store relDat
  piDat=reshape(relDat(:,1:2,:),2*Nsm,[]);
end
piDat=[piDat; squeeze(abs(mean(exp(1i*crdDat(:,3,:)))))']; %adding in orientation variation measure
if(exist('symmStore','var') && ~symmStore); return; end
%% Store computed data
if(exist('stPts','var')); 
  piDatAll=cat(2,piDatAll,piDat(:,stPts));
else; piDatAll=cat(2,piDatAll,piDat(:,stIx:stRes:end));
end
return