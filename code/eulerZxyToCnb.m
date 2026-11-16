function Cnb=eulerZxyToCnb(attiN)
%{
for Cnb,n is up,b is down
given n,we obtain b through rotation with the attitude angles in sequence.
%}

roll  = attiN(1)*180/pi;
pitch = attiN(2)*180/pi;
yaw  = attiN(3)*180/pi;

Rx = rotx(pitch);
Ry = roty(roll);
Rz = rotz(-yaw);

Cnb =Rz*Rx*Ry;
end