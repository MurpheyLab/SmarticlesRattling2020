%return coordinates of 4 points from smarticle position 5-tuple
%in: [cx,cy, theta, al1, al2] - c.o.m coord, body orientation, arm angles (left, right)
%out: [x1,y1, x2,y2, x3,y3, x4,y4]
function [crd] = smcle2coord(sm) %#codegen
global A B;
delB=B/2*[sin(sm(:,3)), -cos(sm(:,3))];
crd=zeros(length(sm(:,1)),8);
crd(:,3:4)=sm(:,1:2)-delB; crd(:,5:6)=sm(:,1:2)+delB; 
crd(:,1:2)=crd(:,3:4)+ A*[-sin(sm(:,3)-sm(:,4)), cos(sm(:,3)-sm(:,4))];
crd(:,7:8)=crd(:,5:6)- A*[-sin(sm(:,3)+sm(:,5)), cos(sm(:,3)+sm(:,5))];
end