%response when force applied to the right arm or corner
%in: arm angle, force angle (to smcle body), arm length (in units of B)
%out: [px,py,f] = pivot coordinates (in units of B) and critical force
%required to move (in units of friction force; f>0 => CCW rotation)
function [out]=pivotArm(al,th,A) %#codegen
    secth=sec(th); tanth=tan(th);
    %the constant out front quantifies how far off the pivot can be (larger=>further away). "right" value: 0.081
    py=0.081./(-2*A*secth.*sin(al-th)+tanth); %approximation for py
    pylnpy2=-0.6*atan(10*py); %differentiable approximation for py*log(py^2) for py in (-0.5,0.5)
    px=tanth/2.*pylnpy2;
    f=-secth.*pylnpy2;
    if(abs(px)>50); px=sign(px)*50; end
    if(abs(py)>50); py=sign(py)*50; end %cap pivot point to avoid divergences
    out=[px,py,f];
end