%resolve all collisions sequentially and iteratively until clear, returning new sm position values
%in: list of smcle 5-tuples, how much to move sm-cles (0-just check collisions), 
%friction coefficients, order of parallel links, ring crds, ring friction (0-no ring)
%out: T/F if just checking, list of new 3-tuples (cx,cy,theta) if resolving
%parOrd - ordering of parallel links, to be kept for a few ticks; xRg -
%ring coordinates (x,y,r^2)
function [xSm,parOrd, xRg]=resolveCollisions(sm, resolveFl, fricCoeff, parOrd, ring, fricR) %#codegen
global A B tRes; Nsm=length(sm(:,1));
smIn=sm;
collFl=true; whCnt=1; xRg=[0,0]; collLks=zeros(Nsm,3); %collision counter
fricSc=1.05; %riction exponential scaling factor
while(collFl) %sequentially resolve all collisions until none left
  collFl=false;
  iterOrd=randperm(Nsm); %randomize each pass to enhance stability
%         crd=smcle2coord_mex(sm);
%         cla; plot(crd(:,1:2:end)',crd(:,2:2:end)','-','LineWidth',2); axis([-0.5,0.5,-0.5,0.5]*(5+2*B)*1.3); axis square;
  
  for smi=1:Nsm
    smi=iterOrd(smi);
    %----Interaction with ring-------
    if(fricR>0) %determines if ring is used
      [dSm,dRg,collTmp]=pushBoundary(sm(smi,:)-[ring(1:2),0,0,0],ring(3),fricCoeff(smi),fricR);
      collLks(smi,:)=collLks(smi,:)+collTmp; %collisions with boundary
      ring(1:2)=ring(1:2)+dRg; sm(smi,1:3)=sm(smi,1:3)+dSm;
      if(any(dRg)); fricR=fricR*fricSc; end %to stop the ring moving eventually
      if(any(dSm) || any(dRg)); collFl=true; 
        fricCoeff(smi)=fricSc*fricCoeff(smi);
      end
    end
    %--------------------------------
    %check for potential collisions based on distance to other smarticles:
    smC = find(sum((repmat(sm(smi,1:2),Nsm,1)-sm(:,1:2)).^2,2) < (2*A+B).^2);
    for ci=1:length(smC)
      smci=smC(ci); %index of the colliding smcle
      if(smci==smi); continue; end %no collision with itself
      %             if(resolveFl) %plot current pair being checked:
      %                 crd=smcle2coord(sm([smi,smci],:));
      %                 cla; plot(crd(:,1:2:end)',crd(:,2:2:end)','-o','LineWidth',2); axis([-0.5,0.5,-0.5,0.5]*12);
      % %                 pause();
      %             end
      
      [res,parOrd(:,:,smi,smci),flipPar, collTmp]...
        =pairCollision(sm([smi,smci],:),resolveFl,fricCoeff([smi,smci]),parOrd(:,:,smi,smci)); %resolve pair collision
      collLks([smi,smci],:)=collLks([smi,smci],:)+collTmp;
      if(any(any(flipPar))); parOrd(:,:,smci,smi)=(flipPar.*parOrd(:,:,smi,smci))'; end
      if all(all(~res)); continue; end %if no collision, continue
      if whCnt>100*Nsm; 'Warning: cannot resolve collision' 
        continue; end %if can't resolve, just move on
      collFl=true; %else set collision flag
      if(all(sm(smi,1:3)==res(1,:))) %if mv==2
        fricCoeff(smci)=fricSc*fricCoeff(smci); %then increase the friction of 2 s.t. it's unlikely to move again
      else;  fricCoeff(smi)=fricSc*fricCoeff(smi);
      end
      sm([smi,smci],1:3)=res;
      sm12v=sm(smi,1:2)-sm(smci,1:2); n12v=norm(sm12v); %vector between smarticles
      if(n12v<0.1*B || (whCnt>10*Nsm && mod(whCnt,3)<1)) %can't resolve collisions
        sm([smi,smci],1:2)=sm([smi,smci],1:2)+...
          [sm12v;-sm12v]/n12v*rand*0.01*B./[fricCoeff([smi,smci]),fricCoeff([smi,smci])]; %rotate COM; add repulsive noise to prevent sticking together
      end
    end
  end
  if(~resolveFl); xSm=collFl*ones(Nsm,3); return; end %if just checking collisions
  %       collFl=false;
  whCnt=whCnt+1;
  if(max(fricCoeff)>50 || any(max(abs(smIn(:,1:3)-sm(:,1:3)))>[B,B,1]*tRes*6*resolveFl)) %too hard to resolve - jam a motor
    xSm=collLks; return;%return collision count instead of new coordinates 
  end
end
%     whCnt
%     fricCoeff
xSm=sm(:,1:3);
xRg=ring(1:2);
end