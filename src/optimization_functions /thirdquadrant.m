function [an]=thirdquadrant(angle)
angle = mod(angle, 360);

if angle>=0 && angle<90 %1->3
    an=angle+180;

elseif angle>=90 && angle<180 %2->3
    an=360-angle;

elseif angle>=180 && angle<270 %3->3 
    an=angle;

elseif angle>=270 && angle<360 %4->3
    an=540-angle;
end
end