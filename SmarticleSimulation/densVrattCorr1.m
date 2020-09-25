% Pavel Chvykov
%Analyzing smarticle data, plotting local density at s.s., vs local
%rattling value to look for a correlation
%also plotting raw data as scatter plot in projected down configuration
%space
%clean up data:-----------
% wavDatAll(isnan(wavDatAll))=1E-10;
% piDatAll(isnan(piDatAll))=1E-10;
% crdDatAll(isnan(crdDatAll))=1E-10;
% wavDatAll(abs(wavDatAll)>1E100)=1E-11;
% piDatAll(abs(piDatAll)>1E2)=1E-11;

oEx=4; nEx=4;%[4,14]-1; %nEx=oEx;
crdCorr=round(2*prd/stRes); %time-segment for R calc
projSS=false(1,size(piDatAll,1)); projSS(1,:)=true; nDim=sum(projSS); %subspace projection %[1:5,Nsm+(1:5)]%[1:5,Nsm+(1:5)]
%% Pull out the data arrays
rattDensSt=[]; datIx=[]; ordDensSt=[]; ordDatSt=[]; corrEigSt=[]; figure;
% =======Steady-state data: ordDat========
ordDat=piDatAll(projSS,(tAll(:,3))==oEx); ordDatAll=ordDat; tDat=tAll((tAll(:,3))==oEx,1); 
selOP=false(1,size(ordDat,2)); %U-shaped configurations
%in regular intervals:
oSt=(3*prd/stRes);  selOP(round(1:oSt:size(ordDat,2)-oSt))=true; % ordered points sampling 
% selOP(oSt-50+1:oSt:end-oSt)=true; %if using many medium runs
%according to arms:
% oSt=35; selOP=squeeze(mean(crdDatAll(:,4:5,tAll(:,3)==oEx))); selOP=(all(abs(selOP-[1;1]*pi/2)<0.3)); selOP(end-100:end)=false; is=1; %for random motion - select in-phase points
% while is<=length(selOP); if(selOP(is)); selOP(is+1:is+oSt-1)=false; is=is+oSt; else; is=is+1; end; end

tmp=crdDatAll(:,4:5,tAll(:,3)==oEx); clf; plot(reshape(tmp(:,:,selOP),Nsm*2,[])'); ylim([-2,2]); drawnow; %to ensure selOP is ok
selOPss = selOP & tDat'>0.3.*max(tDat); %cut off transients
nullDat=ordDat(:,selOP); %will overwrite if shortWavs
ordDat=ordDat(:,selOPss);

% =======Full sampling: nullDat===========
nCyc=3; nRoll=1; tSh=0; %cycles per roll-out, number of roll-outs 
%   ppp=diff(find(diff(tAll(:,2)) & floor(tAll(2:end,3))==nEx,2)); nCyc=round(ppp/prd*stRes); %single run length (for short runs)
runLen=round(nCyc*prd/stRes); ppp=nRoll*runLen;
nullDat=piDatAll(projSS,any(tAll(:,3)==nEx,2));
nullDatFull=reshape(nullDat,nDim,runLen,nRoll,[]);
nMx=size(nullDat,2); selIx=tSh+1:ppp:nMx; 
nullDat=nullDat(:,selIx); nPts=length(nullDat(1,:));
nullDIx=(tAll(any(tAll(:,3)==nEx,2),:)); nullDt=(nullDIx(selIx,1)); nullDIx=(nullDIx(selIx,3)); %index showing which dIx the point came from

% =======Seed points: seedDat=============
sPts=min(1000,size(nullDat,2)); seedSel=randsample(size(nullDat,2),sPts); seedDat=nullDat(:,seedSel); %subsample
% cntPtsDat=[cntPtsDat,nullDat(:,randsample(size(nullDat,2),600))];continue;
seedDIx=nullDIx(seedSel); seedDt=nullDt(seedSel);


%% Calc Rattling for each nullDat
velAll=(nullDatFull-nullDatFull(:,1,:,:))./sqrt(0:runLen-1); velAll=velAll(:,3:crdCorr,1:end,:); %since diffusion ~sqrt(t)
velAll=permute(reshape(velAll,nDim,[],size(velAll,4)),[1,4,5,3,2]); %nDim x 1 x 1 x nPts x crdCorr
velAll=velAll-mean(velAll,5); %center velocoties (is this necessary?)
% stepDist=squeeze(mean(sum(velAll.^2,1),5)); %find mean velocity 
corrEig=zeros(nDim,nPts); corrSt=zeros(nPts,nDim,nDim);
w=waitbar(0); 
for ti=1:nPts %compute rattling of each point
  corr1=mean(velAll(:,:,:,ti,:).*permute(velAll(:,:,:,ti,:),[2,1,3,4,5]),5);
  corrEig(:,ti)=abs(eig(corr1)); corrSt(ti,:,:)=corr1;
  
  if(mod(ti,nPts/50)<1); waitbar(ti/nPts,w,'ratt'); end
end
close(w);
%% Average density and rattling

aveVol=zeros(1,sPts); nullDens=aveVol; ordDens=aveVol; nOrdPts=aveVol; rattDens=aveVol;
nnn=ceil(sum(projSS)*3); %if(sum(projSS)<=10); nnn=15; else; nnn=50; end
w=waitbar(0);
for ip=1:sPts 
  %===== Density========
  currPt=seedDat(:,ip); %center of the cell being examined
  tmpDiff=(nullDat-currPt); nullDist=sum(tmpDiff.^2);
  [nullDist, nrst]=sort(nullDist); %find nnn nearest points 
  tmpDiff=tmpDiff(:,nrst(1:nnn)); %tmpDiff=tmpDiff-mean(tmpDiff,2); currPt=mean(nullDat(:,nrst(1:nnn)),2); %shift to cloud center
%   plot(nrst(1:nnn)); ylim([0,length(nullDist)]); pause(0.01); %check if all nrst nbr points are sequential in time
  corrMxLoc=tmpDiff*tmpDiff'/nnn; aveVol(ip)=sqrt(det(corrMxLoc)); %volume that nn points take up - sqrt det of corr matrix
  nullDens(ip)=sum(nullDIx(nrst(1:nnn))==nEx(1))./aveVol(ip);%./nullDist(nnn)^(sum(projSS)/2); %null density, for normalization
  tmpDiff=(ordDat-currPt); 
  ordDist=(sum(tmpDiff.*(corrMxLoc\tmpDiff))); nOrdPts(ip)=sum(ordDist<5)./length(ordDist); %fraction of ordDat in the nbhd
  ordDens(ip)=nOrdPts(ip)./aveVol(ip); %density

  %===== Rattling ========
%   corrAve=squeeze(mean(corrSt(nrst(1:nnn),:,:),1)); rattDens(:,ip)=eig(corrAve); %average T_eff matrices (bad)
  rattDens(ip)=mean(log(prod(2*pi*exp(1)*corrEig(:,nrst(1:nnn)))))/2; %average rattling values
  
  if(mod(ip,sPts/50)<1); waitbar(ip/sPts,w,'pdf'); end
end
close(w);

rattDensSt=[rattDensSt,rattDens]; datIx=[datIx,oEx.*ones(1,length(rattDens))]; ordDensSt=[ordDensSt,ordDens];
ordDatSt=[ordDatSt,ordDat]; corrEigSt=[corrEigSt,corrEig];

%% piDat Scatter plots
figure; ax1=subplot(122);
xFit=log(prod(2*pi*exp(1)*corrEig(:,:)))/2; xFit=rattDens; %xFit=log(prod(2*pi*exp(1).*rattDens(:,:)))/2;
% xFit=log(ordDens);
if(sum(projSS)>300); plD=[1,2,3; 4,5,6; 7,7,7]';plD=[2,1,3];%plD=[1,4,7]; 
else; plD=[1,2,3]; end %plD=[6,5,4]; nullDat=eV\nullDatB; ordDat=eV\ordDatB;
scatter3(-mean(seedDat(plD(:,1),:),1),-mean(seedDat(plD(:,2),:),1),mean(seedDat(plD(:,3),:),1)...
  ,0.1,xFit(:))%,'MarkerEdgeAlpha',.4)
axis([-1.5,1,-1.5,1,-1,1.5]); view([58,31]); colormap jet;
axis([-1.5,1,-1.5,1,-1,1.5]); view([-149-90,20.5]);
colormap jet; %caxis([-24,-16]);%axis([-1,1,-1,1,0,0.8]); view([-1/2,-1,1/2]);
set(gca,'xticklabels',[],'yticklabels',[], 'zticklabels',[]); %axis vis3d; view([0.5,0.5,0.5]); %colorbar;
% xlabel('x_1','FontSize',16); ylabel('x_2','FontSize',16); zlabel('x_3','FontSize',16); 

ax2=subplot(121); %scatter3(nullDat(plD(1),:),nullDat(plD(2),:),nullDat(plD(3),:),1,log(ordDens1));
scatter3(-mean(ordDat(plD(:,1),:),1),-mean(ordDat(plD(:,2),:),1),mean(ordDat(plD(:,3),:),1),0.1);%[0,0,0]); 
% hold on; scatter3(nullDat(plD(1),:),nullDat(plD(2),:),nullDat(plD(3),:),ptSz)%,log((rattDensB)))%,'MarkerEdgeAlpha',.4)
set(gca,'xticklabels',[],'yticklabels',[], 'zticklabels',[]); %axis vis3d;
linkprop([ax1,ax2],{'CameraPosition','CameraUpVector'}); rotate3d on
axis([-1.5,1,-1.5,1,-1,1.5]); view([58,31]); colormap jet;
axis([-1.5,1,-1.5,1,-1,1.5]); view([-149-90,20.5]);
% xlabel('x_1','FontSize',16); ylabel('x_2','FontSize',16); zlabel('x_3','FontSize',16); 
% view([0.5,0.5,0.5]); while true; camorbit(ax1,1.5,-0.); camorbit(ax2,1.5,-0.); drawnow; end %rotate 3D
% hold on; ptIx=8; plot3(-mean(nullDat(plD(:,1),ptIx),1),-mean(nullDat(plD(:,2),ptIx),1),mean(nullDat(plD(:,3),ptIx),1),'*k');
% return
%% p_ss v det
ordDens1=ordDens; ordDens1(ordDens==0)=min(ordDens(ordDens>0));
xFit=log(prod(2*pi*exp(1)*corrEig(:,seedSel)))/2; xFit=(rattDens)*1; %xFit=log(prod(2*pi*exp(1)*rattDens(:,:)))/2; nullDIx(end)=0;
figure; SCAT=scatter(xFit,log(ordDens1),1.,seedDIx);%(nullDens./maxplko(nullDens)*1)); %colormap(flip(gray));
colormap(jetR); %axis equal; %axis([-27,-6,4,28]); %axis([-74,-49,4,24]);
% xlabel('$\mathcal{R}$','FontSize',16,'Interpreter','latex');%,'FontWeight','bold'); 'log(det\langleT_{eff}\rangle)'
% ylabel('log(p_{s.s.})','FontSize',16);%,'FontWeight','bold');
% hold on; plot(xFit(ptIx),log(ordDens1(ptIx)),'*k')
yFit=log(ordDens); keep=yFit>-Inf & xFit>-80; xFit=xFit(keep); yFit=yFit(keep); tmp=corrcoef(xFit,yFit); tmp1=polyfit(xFit,yFit,1); title(['slope=', num2str(round(tmp1(1),2)),',  \rho=',num2str(-round(tmp(1,2),2))],'FontSize',16);
return
%% p_ss v trace
xFit=log(sum(2*pi*exp(1)*corrEig(:,seedSel)));
figure; SCAT=scatter(xFit,log(ordDens),1.,seedDIx);%(nullDens./max(nullDens)*1)); %colormap(flip(gray));
colormap(jetR); %axis equal;
xlabel('|v^2|','FontSize',16,'FontWeight','bold'); %log(tr\langleT_{eff}\rangle)
ylabel('log(p_{s.s.})','FontSize',16,'FontWeight','bold');
yFit=log(ordDens); xFit=xFit(yFit>-Inf); yFit=yFit(yFit>-Inf); tmp=corrcoef(xFit,yFit); tmp1=polyfit(xFit,yFit,1); title(['slope=', num2str(round(tmp1(1),2)),',  \rho=',num2str(-round(tmp(1,2),2))],'FontSize',16);
return
