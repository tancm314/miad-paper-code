function atti=cnbToEulerZxy(Cnb)

%{
revised by Tan Caiming on 15 April 2021 in NJUPT
make it robust
Cnb=rotz(-yaw)*rotx(pitch)*roty(roll)
%}
Cbn=Cnb';
if(abs(Cbn(2,3))>1)
    sp=sign(Cbn(2,3));
    disp('rotatin Matrix is not formal!')
else
    sp=Cbn(2,3);
end
pitch=asin(sp);        %pitch -90~90
cp=cos(pitch);

if(abs(cp)<1e-20)
    disp('roll and yaw cannot be determined when pitch is +/-90°!');
end
roll=atan2(-Cbn(1,3),Cbn(3,3));  %roll -180°~180°
yaw=atan2(Cbn(2,1),Cbn(2,2));   %yaw 0~360°
if(yaw<0)
    yaw=2*pi+yaw;
end
atti=[roll;pitch;yaw];
