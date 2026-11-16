function Rv=quaternionToRotationVector(Q)
%{
Revised by Caiming Tan on 23 March 2023
%}
a=2*atan2(norm(Q(2:4)),Q(1));
% n=Q(2:4)/norm(Q(2:4));
if sin(a/2)~=0
    n=Q(2:4)/sin(a/2);
else
    n=[0;0;0];
end
Rv=a*n;