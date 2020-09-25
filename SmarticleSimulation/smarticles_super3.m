% Smarticles simulation
% Pavel Chvykov
% clear all;
%-------------------------
rng(5001);
global A B tRes latFric; B=1; A=0.9; %smarticle size %B1; A0.9
latFric=0.2; %set at 0.2 (increasing makes things less stable)
Nsm=3; tRes=2*pi/200; %0.03; %number of smarticles
prd=2*pi/tRes; stRes=round(prd/50);
nCyc=3*10; runLen=round(nCyc*2*pi/tRes); %length(t); %length of one run %3000 for long runs
t=(0:1*runLen)*tRes; nRuns=50; %0:tRes:6*pi;%15000; %8*pi %time domain - repeated runs
%Clear / initialize data storage arrays:
crdDatAll=[]; piDatAll=[]; tAll=[]; fnameAll=[]; ringAll=[]; 

windSize=(B+2*A)/3*sqrt(Nsm)*1.05; %1.05 - for N=3 %%exp't: 19.2/5.2*6/7=3.165
scRat=[20,1];  
fricR=1.3E6; dragRv=0*B/500; %ring friction (fricR=0 -> no ring, confining potential instead)
inertCoeff=0*exp(-tRes/0.05); resDist=100; %inertia: 0:none, 1 :full; max sm-cle displacesment per t-step (~1) %exp(-tRes/0.01);
livePlot=1.; plFrom=0.*pi; fpp=24; %figure; pause %plot flag, lines thickness, frames per period
plScl=1; plRange=2; %live plot: line thickness and display range
T=0.E-3; armAmpRnd=0.01; armSpdRnd=0.0; %additional noise sources
% paramAll=[2:15,2:15,2:15];%0.1:0.1:1;%exp(0:0.3:3.3);%inertial coefficient
% load('191031nullSample.mat'); %null sampling of config space
% nullCrd=randsample(size(nullSample,3),nRuns); %sample null distribution
for ordIC=false%[true,false] %for short-wave runs - check nOrdRuns first
  if(ordIC)
  oIx=14; ordCrd=(crdDatAll(:,:,tAll(:,3)==oIx)); %ordered steady-state samples (oIx - for long runs)
  sel=squeeze(mean(ordCrd(:,4:5,:))); sel=(all(abs(sel-pi/2)<0.3)); is=1; ppp1=1*round(prd/stRes);%for random motion - select in-phase points by arm angles
  while is<=length(sel); if(sel(is)); sel(is+1:is+ppp1-1)=false; is=is+ppp1; else; is=is+1; end; end
  ordCrd=ordCrd(:,:,sel); 
%   ordCrd=ordCrd(:,:,1:prd/stRes:end); 
%   ptsUsd=randsample(size(ordCrd,3),nRuns);
  ordCrd=ordCrd(:,:,unique(ceil((1:nRuns)/nRuns*size(ordCrd,3)))); nRuns=size(ordCrd,3);
  size(ordCrd,3)
  pause;
%   ordCrd=cat(3,crdDatMax(:,:,1:round(nRuns/4)),crdDatMax(:,:,end-round(nRuns/4):end));
%   ordCrd=cat(3,ordCrd,crdDatMax(:,:,randsample(round(nRuns/4):size(crdDatMax,3)-round(nRuns/4),nRuns-size(ordCrd,3))));
  end
%   gatesListAll=[];
for expIx=[4] %experiment index: determines gait and labels stored data, as follows
%expIx selects different gaits: 1=random gait; 2=correlated random gait; 3=square gait; 
%4=distinct periodic gait of period 3 (set seed); expIx > 100: distinct periodic gait of period 4
%expIx in range (4,5]-add drive randomness; 20=l3gt33; 21=mixGt(l3gt33+seed); in (20,21)-add inertia to 20

%   nCyc=round(3*60/4*paramAll(expIx)); runLen=round(nCyc*2*pi/tRes); t=(0:1*runLen)*tRes; %for changing gait period
  if(ordIC); expIx=expIx+0.1; end %to distinguish the label

if(expIx>20 && expIx<21); inertCoeff=exp(-tRes./(expIx-20)); end%exp(-tRes./paramAll(expIx));
for repIx=1:nRuns%2000 %10 %loop to just repeat the runs
  [expIx, repIx, ordIC]
%   rng(16);
%======Set motion gaits=================
fricCoeff=1E-1*(1+0.*(rand(Nsm,1)-0.5)); %friction coefficients for different smcles

Nwall=0;
%% Define the Gait:----------------------
rng(33);%rng(8);%rng(27); %seed for the choice of random gait
tStep=round(2*pi/4/tRes); maxV=1.8; nSeg=4E4; %Random gait %int: (randi(2,Nsm,2,1,nSeg)*2-3)
halfFl=false; 
if(expIx>=4);nSeg=3;end
if(expIx>=100);nSeg=4;end
while(~halfFl) 
% nSeg=paramAll(expIx);  
if(expIx<=1 || expIx>=4)
patt=(randi(2,Nsm,2,1,nSeg)*2-3).*ones(Nsm,2); patt(:,:,1,nSeg:4*nCyc:end)=1; %unique motion patterns for each Smcle
nRet=12; if(nSeg>nRet);patt(:,:,1,nRet:nRet:end)=1; end %return back to u-shape to allow fair comparison
elseif(expIx==2); rng('shuffle');
  patt=(randi(2,1,2,1,nSeg)*2-3).*ones(Nsm,2); patt(:,:,1,nSeg:4*nCyc:end)=1; %patt(:,:,1,24:24:end)=1; %unique motion patterns for each Smcle
nRet=12; if(nSeg>nRet);patt(:,:,1,nRet:nRet:end)=1; end %return back to u-shape to allow fair comparison
elseif(floor(expIx)==3); patt=(repmat(permute([1,-1; -1,-1; -1,1; 1,1],[4,2,3,1]),[Nsm,1,1,1])); nSeg=size(patt,4); %square gait
end
moves=diff(cat(4,patt,patt(:,:,1,1)),1,4); moves=sum(moves(:)~=0)/numel(moves); halfFl=0.45<moves && moves<0.65; %ensure moving 50% of the time
end
%for experimental code:
pattExp=reshape(permute(patt.*[-1,1],[2,1,4,3]),[],nSeg).*90+90; pattExp=[pattExp(:,end),pattExp(:,1:end-1)];

  % patt=(repmat(permute([-1,-1; -1,1; 1,1],[4,2,3,1]),[Nsm,1,1,1])); %triangle gait
% patt=repmat([1,1],[Nsm,1,1,nSeg]); %hold constant position
% nSeg=8; patt=ones(Nsm,2,1,nSeg); patt(1,1,1,1:3)=-1; patt(2,1,1,2:4)=-1; patt(3,1,1,3:5)=-1; patt(:,2,:,:)=patt([2,3,1],1,:,:); %staggered arm motion
% gatesList=repmat(reshape(patt.*ones(1,1,tStep),Nsm,2,[]),[1,1,ceil(length(t)./tStep/nSeg)]); gatesList(:,:,1)=ones(Nsm,2); gatesList=pi/2*gatesList;
patt33=[]; patt33(:,:,1,3)=ones(3,2); patt33(:,:,1,1)=[-1,-1; -1,1; -1,-1]; patt33(:,:,1,2)=[-1,-1; 1,1; 1,-1];
patt34=[]; patt34(:,:,1,3)=ones(3,2); patt34(:,:,1,1)=[-1,1; 1,-1; -1,1]; patt34(:,:,1,2)=[-1,-1; -1,-1; -1,1];
pattNW1=[]; pattNW1(:,:,1,3)=ones(3,2); pattNW1(:,:,1,1)=[-1,-1; 1,1; 1,-1]; pattNW1(:,:,1,2)=[-1,-1;-1,1; -1,-1];
if(expIx==14); %switching gait
  gatesList=cat(4,repmat(patt33,[1,1,1,40]),repmat(patt34,[1,1,1,30]),...
                  repmat(patt33,[1,1,1,15]),repmat(patt34,[1,1,1,45]),...
                  repmat(patt33,[1,1,1,20]),repmat(patt34,[1,1,1,33]));
elseif(expIx==20); gatesList=repmat(patt33,[1,1,1,ceil(length(t)./tStep/nSeg)]);
elseif(expIx==21); gatesList=[]; pIx=0; %index which gait to run
  for ip=1:ceil(length(t)./tStep/nSeg)
    if(pIx); gatesList=cat(4,gatesList,patt33); 
    else; gatesList=cat(4,gatesList,patt);
    end
    if(rand<0.2); pIx=abs(pIx-1); end %sometimes switch gait
  end
else; gatesList=repmat(patt,[1,1,1,ceil(length(t)./tStep/nSeg)]);
end

if(expIx>4 && expIx<=5)
% rndChg=randsample(numel(gatesList),round((0.1)*numel(gatesList))); gatesList(rndChg)=-gatesList(rndChg); gatesList(:,:,:,3*nSeg:3*nSeg:end)=1; %randomly flip a fraction of the moves %10.^(2.5*(expIx-4))
rndChg=randsample(numel(gatesList),round(10.^(2.*(expIx-5))*numel(gatesList))); gatesList(rndChg)=2*randi(2,numel(rndChg),1)-3;  %randomly change some of the moves
end
gatesList(:,:,:,4*nSeg:4*nSeg:end)=1; %to allow for fair comparison across gaits
gatesList=reshape(gatesList.*ones(1,1,tStep),Nsm,2,[]); 
gatesList(:,:,1)=ones(Nsm,2); %start all smarticles in U-shape
gatesList=(pi/2-0*0.25)*gatesList; %one move = pi/2 in time
scl=reshape(repmat(permute(1.0-armAmpRnd*(rand(Nsm,2,ceil(size(gatesList,3)/tStep/2))-0.5),[1,2,4,3]),[1,1,2*tStep,1]),Nsm,2,[],1);
gatesList=gatesList.*min(1,scl(:,:,1:size(gatesList,3)));

% figure; plot(squeeze(gatesList(1,1,:)),squeeze(gatesList(1,2,:)))
% figure;plot(squeeze(gatesList(:,1,1:20*pi/tRes))')
% return
%%
rng(repIx*10+1E6*expIx+122); %change IC for new runs
xSm=zeros(Nsm,5);
xSm(:,4:5)=gatesList(:,:,1); %gates(0);
%note: fric grows exp with each collision resolution step and fric>50 jams motor
% fricCoeff=[1.5,0.4];
%--------Set wall--------------
% xSm(1:Nwall,1:3)=zeros(Nwall,3);
% xSm(1:Nwall,2)=((1:Nwall)-3)*(B+2*A)*0.6;
% xSm(1:Nwall,3)=pi/2;
% out=1.5;shft=0.05; xSm(1:Nwall,1:3)=[shft,out,pi/2; out,-shft,0; -shft,-out,-pi/2; -out,shft,pi];
% xSm(1:Nwall,4:5)=-ones(Nwall,2)*pi/6;
% fricCoeff(1:Nwall)=1E10;

%% ======Set smarticle initial positions: [cx,cy, theta, al1, al2] - c.o.m coord, body orientation, arm angles
% xSm=[2,1,pi/4,pi/8,pi/8];
% xSm(:,1:3)=(rand(Nsm,3)-0.5).*repmat([10,10, 2*pi],Nsm,1);
if(ordIC) %sample from p_ss
  xSm=ordCrd(:,:,repIx);
  xSmIC=xSm; %store IC to repeat runs
elseif(exist('nullSample','var')) %use null sampling from file
  xSm=nullSample(:,:,nullCrd(repIx));
else %randomly drop smarticles and resolve collisions
for smi=Nwall+1:Nsm 
    collFl=true;
    ii=1;
    while(collFl)
        ii=ii+1;
        collFl=false;
        xSm(smi,1:3)=(rand(1,3)-0.5).*[windSize,windSize, 2*pi]; %add new smcle
        for smic=1:smi-1 %check that it doesn't collide with previous ones
            [res,~,~]=pairCollision(xSm([smi,smic],:), false, [1,1], zeros(9,1));
            collFl=collFl || res(1); %only check for collisions
        end
        if(ii>10000); error('cannot build IC'); end
    end
end
end

% Set specific initializations:--------------
% xSm(:,1:3)=[-3,0,pi/2; -1.1,0.1,pi/2-0.1];
% xSm(:,1:3)=[-0.8,0,0; 0,0.1,pi; 0.8,0,0]; %+-+ IC state
% xSm(:,1:3)=[-0.6,0,0; 0.3,0.5,-2*pi/3; 0.3,-0.5,2*pi/3]+[0,0,pi]; %rot-symmetric IC
% xSm(:,1:3)=[0,2.,pi/2; 0,0,pi/2; 0,-2.,pi/2]; %stacked IC - demonstrate gait
%%-------- From exp't tracked data---------------------
% iMov=1; crdDat=[];
% crdDat(:,1,:)=permute(movs(iMov).x,[2,3,1])*100/5.2/1.; %in units of smarticle size
% crdDat(:,2,:)=permute(movs(iMov).y,[2,3,1])*100/5.2/1.;
% crdShift=mean(mean(crdDat(:,1:2,:)),3);
% crdDat(:,1,:)=crdDat(:,1,:)-crdShift(1);
% crdDat(:,2,:)=crdDat(:,2,:)-crdShift(2);
% crdDat(:,3,:)=permute(movs(iMov).rot,[2,3,1])+pi/2;
% xSm=crdDat(:,:,1); xSm(:,4:5)=gates(0);
% 
% crd=smcle2coord_mex(xSm);
% clf; plot(crd(:,1:2:end)',crd(:,2:2:end)','-','LineWidth',2); axis([-0.5,0.5,-0.5,0.5]*(windSize+2*B)*2.); axis square;
% return
%% ====Confining potential=======
rRad=windSize; ring=[0,0,rRad^2]; %hard ring bdry coordinates, size
%Soft - Bowl potential
Upow=20;
% Fcent = @(vec) -vec.*min((vec(:,1).^2+vec(:,2).^2).^(Upow/2-1).*Upow/windSize^(Upow-2),...
%     2)/windSize^2;
%Circular "plate" potential
Fcent = @(vec) -0*15*(vec)/B.*subplus(1-rRad./vecnorm(vec,2,2));
%Rectangular plate:
% Fcent = @(vec) -vec.*scRat/B.*subplus(1-windSize./max(abs(vec.*scRat),[],2));

%% ======Simulation Run=============
% hold on;
% nextAl=gates(t(1));
crdDat=zeros(Nsm,5,length(t)); ringDat=zeros(3,length(t));
parOrd=zeros(3,3,Nsm,Nsm); %store ordering of parallel links (to avoid passing through)
%  w=waitbar(0,['run ',num2str(repIx)]);
crdDat(:,:,1)=xSm; ringDat(:,1)=ring; tCol=0; delAl=0; clf;
for ti=2:length(t) %main time ticks loop
%   if(mod(ti,pi/tRes)<1); A=0.8+0.2*rand; end
  lastAl=xSm(Nwall+1:end,4:5); 
  nextAl=gatesList(:,:,ti);%gates(t(ti));
  delAl0=delAl; delAl=(nextAl(Nwall+1:end,:)-lastAl);
  if(maxV>=0); delAl=sign(delAl).*min(abs(delAl),maxV*tRes*ones(Nsm-Nwall,2).*abs(1+armSpdRnd*randn(Nsm,2))); %move to new angle or at max speed
  else; maxV=-5; delAl=(1-2*sqrt(-maxV).*tRes)*delAl0 - tRes.*tRes.*maxV.*(delAl); %spring-like pull to new position, with inertia
  end
  xSm(Nwall+1:end,4:5)=xSm(Nwall+1:end,4:5)+delAl; %move the arms
%   ring(3)=rRad2*(1+0.2*rand);  
  ring(1)=ring(1)+dragRv*tRes; %pull the ring
    while(true)
      if(ti==2); [res,parOrd,ring(1:2)]=resolveCollisions(xSm,1000, fricCoeff, parOrd, ring, fricR); %resolve all collisions sequentially
%         rng('shuffle'); res=res+0.1*(rand(Nsm,3)-0.5); %perturb IC a bit
      else
%       if(t(ti)<583) %for Debugging
        [res,parOrd,ring(1:2)]=resolveCollisions_mex(xSm,resDist, fricCoeff, parOrd, ring, fricR); %resolve all collisions sequentially
%       else; [res,parOrd,ring(1:2)]=resolveCollisions(xSm,true, fricCoeff, parOrd, ring, fricR); %resolve all collisions sequentially
%       end
      end
      if (all(all(mod(res,1)==0))) %couldn't resolve collision - jam motor & try again
        res=[res(:,1),res(:,3)]; res(abs(xSm(:,4:5)-lastAl)<1E-5)=0; [mxColl,smIx]=max(res); [~,armIx]=max(mxColl); smIx=smIx(armIx);
        xSm(smIx,3+armIx)=lastAl(smIx,armIx); %revert the most colliding arm
        if(sum(sum(res>0))<=1); 'Warning: all motors jammed'
          break; end
      else %apply the collision-resolving move
        xSm(:,1:3)=res; break;
      end
    end
        parOrd=parOrd/2; parOrd(abs(parOrd)<0.02)=0; %decay memory of parallel over time
%     xSm(:,1:2)=xSm(:,1:2)+tRes*Fcent(xSm(:,1:2))./fricCoeff;

%     if (~fricR) %confining potential (only if no ring)
      crd=smcle2coord_mex(xSm);
      Ftmp=Fcent(reshape(crd',2,[])'-ring(1:2));
      xSm(:,1:2)=xSm(:,1:2)+tRes*(Ftmp(1:4:end,:)+Ftmp(2:4:end,:)+Ftmp(3:4:end,:)+Ftmp(4:4:end,:))./4./fricCoeff;
%     end
    
    nzAmp=sqrt(tRes./fricCoeff.*[B,B,1]*T); xSm(:,1:3)=xSm(:,1:3)+nzAmp.*randn(Nsm,3); %add randomness
    if(livePlot && mod(ti,2*pi/tRes/fpp)<1 && t(ti)>plFrom) %show smarticles motion live
        crd=smcle2coord_mex(xSm);%pause
        if(livePlot==2);  %overlay - make rainbow plots
          cmap=colormap(hsv_drk); tCol=tCol+1; tmpCol=cmap(ceil(mod(tCol/fpp,1)*63+1E-3),:);
%           plScl=4; if(t(ti)<60); clf; transp=0.2; else; transp=0.05; end %for movie
          hold on; plot(crd(:,1:2:end)',crd(:,2:2:end)','-','LineWidth',2*plScl,'Color',[tmpCol,0.03]); axis([-0.5,0.5,-0.5,0.5]*(windSize+2*B)*plRange); axis square;
%           plot(xSm(:,1)+0.05*cos(xSm(:,3)),xSm(:,2)+0.05*sin(xSm(:,3)),'.','MarkerSize',15,'Color',tmpCol); %can't do transparency here
%           scatter(xSm(:,1)+0.05*cos(xSm(:,3)),xSm(:,2)+0.05*sin(xSm(:,3)),25,tmpCol,'filled','MarkerFaceAlpha',0.3);
        elseif(livePlot==2.01)
          cmap=hsv; tCol=tCol+1; tmpCol=cmap(ceil(mod(tCol/fpp,1)*63+1E-3),:);
          hold on; plot(crd(:,3:2:5)',crd(:,4:2:6)','-','LineWidth',2*plScl,'Color',[tmpCol,0.03]); axis([-0.5,0.5,-0.5,0.5]*(windSize+2*B)*plRange); axis square;
        else %normal plot
        clf; plot(crd(:,1:2:end)',crd(:,2:2:end)','-','LineWidth',2*plScl); axis([-0.5,0.5,-0.5,0.5]*(windSize+2*B)*plRange); axis square;
        hold on; plot(xSm(:,1)+0.05*cos(xSm(:,3)),xSm(:,2)+0.05*sin(xSm(:,3)),'k.','MarkerSize',15*plScl);
        crdShX=crd(:,[3,5]); crdShY=crd(:,[4,6]); jam=abs(xSm(:,4:5)-nextAl)>1E-5;
%         hold on; plot(crdShX(jam),crdShY(jam),'.r','MarkerSize',20); %show jammed motors
%         rush=abs(abs(xSm(:,4:5)-lastAl)-4*tRes)<1E-5; plot(crdShX(rush),crdShY(rush),'.g','MarkerSize',20);
        end
        if(fricR) %plot ring
          th=0:0.1:2.1*pi; plot(windSize*cos(th)+ring(1),windSize*sin(th)+ring(2),'k','LineWidth',2*plScl);
        end
        title([repIx,t(ti)]);
%         if(ti==2); pause; end
        
%         tmp=toc; pause(0.2-tmp); tic; %ensure equal times between frames for video recording
        drawnow;
%         if(mod(tCol,50)==0);pause(0.01);end
        if(livePlot==2.1); return; end
    end
    
    crdDat(:,:,ti)=xSm; ringDat(:,ti)=ring;
    if(~ordIC && ti==2); xSmIC=xSm; end %store IC after resolving collisions
    xSm(:,1:3)=xSm(:,1:3)+(crdDat(:,1:3,ti)-crdDat(:,1:3,ti-1)).*inertCoeff; %Inertia!
%     if(mod(ti,length(t)/30)<=1); waitbar(ti/length(t),w); end
    if(mod(ti,runLen)==2); xSm=xSmIC; %reset to IC
%       nzAmp=sqrt(tRes./fricCoeff.*[B,B,1]*0.2E-4); xSm(:,1:3)=xSm(:,1:3)+nzAmp.*randn(Nsm,3); %perturb IC
      crdDat(:,:,ti)=xSm;
    end

    if(all(any(abs(xSm(:,1:2)-ring(1:2))>rRad,2))) %error - quit
%       crdDat=crdDat(:,:,1:ti); ringDat=ringDat(:,1:ti); t=t(1:ti);%clip time-series
      break; 
    end
end
% close(w); %waitbar
% return

if(ti<length(t)); warning('jumped out of ring: moving to next run');
  continue; end %if quit with error, go on to next run
%% ===========Data Analysis=====================
crdDat(:,:,1)=crdDat(:,:,3); %ensure that collisions resolved in IC
% ppp=20; tSparse=round(1:prd/ppp:length(t)); %Sparse time to ppp - points per period
symmetrize_smarticles1; %compute and store rotation and permutation invariant observables

tAll=[tAll; [t(2:stRes:end)',[repIx,expIx].*ones(size(t(2:stRes:end)'))]];
crdDatAll=cat(3,crdDatAll,crdDat(:,:,2:stRes:end));
% ringAll=[ringAll; ringDat(:,2:stRes:end)'];
end
end
end
%% =========Run data alaysis============
densVrattCorr1; 

return
%% Save generated data
save test1.mat crdDatAll piDatAll tAll tRes A B Nsm latFric prd stRes nCyc runLen nRuns windSize fricCoeff dragRv T armAmpRnd armSpdRnd freqList phaseList Nwall

%% Re-compile _mex files (requires Matlab coder package)
codegen smcle2coord -args {coder.typeof(0,[1000,5],[1,0])}
codegen resolveCollisions -args {coder.typeof(0,[1000,5],[1,0]),1,coder.typeof(0,[1000,1],[1,0]),coder.typeof(0,[3,3,1000,1000],[0,1,1]),[0,0,0],0}
