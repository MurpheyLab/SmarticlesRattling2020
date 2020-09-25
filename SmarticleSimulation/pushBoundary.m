%smarticle interaction with a ring boundary
%in: smcles array (relative to ring center), ring radius^2, smcles frictions, ring friction
%out: dSm=shift of smarticle coordinates, dRg=shift of ring coorinates
function [dSm,dRg,collLks]=pushBoundary(sm,rad2,fricCoeff,fricW) %#codegen
global tRes B A;
crd=smcle2coord(sm); Nsm=length(sm(:,1));
pivAngle=tRes; parShift=tRes*B/3;
dSm=zeros(Nsm,3); dRg=[0,0]; collLks=zeros(Nsm,3);
for si=1:Nsm %loop over smarticles
  outVx=find(crd(si,1:2:end).^2+crd(si,2:2:end).^2>rad2); %vertices sticking out
  collLks(si,1)=collLks(si,1)+sum(outVx==1);
  collLks(si,3)=collLks(si,3)+sum(outVx==4); %store number of collisions on lk
  switch(length(outVx))
    case 0; continue; %no contact with bdry
    case 1 %one contact point - rotate
      [thVx,~]=cart2pol(crd(si,outVx*2-1),crd(si,outVx*2));
      piv=[0,0,0];
      switch(outVx)
        case 1  %pushing with left arm tip
          piv=pivotArm(sm(si,4), -thVx+sm(si,3)-pi/2, A/B); piv([1,3])=-piv([1,3]);
        case 2 %pusing with left corner
          piv=pivotArm(0, -thVx+sm(si,3)-pi/2, 0); piv([1,3])=-piv([1,3]);
        case 3 %if pushing with right corner
          piv=pivotArm(0, thVx-sm(si,3)-pi/2, 0);
        case 4  %if pushing with the right arm tip
          piv=pivotArm(sm(si,5), thVx-sm(si,3)-pi/2, A/B);
      end
      pivF=abs(piv(3).*fricCoeff(si));
      if((pivF-fricW)/(pivF+fricW)<6*(rand-0.5)^3) %move smarticle
        cth=cos(sm(si,3)-pi/2); sth=sin(sm(si,3)-pi/2); rotMx=[[cth,-sth];[sth,cth]]; %rotation matrix
        pivCOM=(rotMx*piv(1:2)'*B)'; %pivot to smarticle COM in lab frame
        dSm(si,1:2)=sign(piv(3))*pivAngle*[1,-1].*flip(pivCOM)./norm(pivCOM);%+trEps*abs(pivAngle)*sm12v/norm(sm12v)*rand;
        dSm(si,3)=sign(piv(3))*pivAngle; %rotate angle
      else %move wall
        dRg=dRg+parShift*[cos(thVx),sin(thVx)];
      end
    otherwise %multiple contact points - translate
      shD=[0,0];
      for vi=1:length(outVx)
        Fvc=crd(si,outVx(vi)*2+(-1:0)); Fvc=Fvc/norm(Fvc);
        shD=shD+Fvc*parShift;
      end
      if((fricCoeff(si)-fricW)/(fricCoeff(si)+fricW)<6*(rand-0.5)^3) %move smarticle
        dSm(si,1:2)=-shD;
        dSm(si,3)=randn*pivAngle;
      else %move wall
        dRg=dRg+shD;
      end
  end
end
end