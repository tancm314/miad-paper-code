function theta=computeAttitudeDistance(atti,atti2,index)
% Revised by Caiming Tan on 14 December 2022 in Hongyunjiayuan community
switch index
    case 1
        Cnb=attitudeAnglesToCnb(atti);
        Cbn=Cnb';
        Cnb2=attitudeAnglesToCnb(atti2);
        Cbn2=Cnb2';
    case 2
        Cbn=atti;
        Cbn2=atti2;
end
C=Cbn*Cbn2';
Q=cnbToQuaternion(C);
v=quaternionToRotationVector(Q);
theta=min(norm(v),2*pi-norm(v));