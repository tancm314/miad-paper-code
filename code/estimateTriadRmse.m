function RMSE_phi=estimateTriadRmse(w1,w2,std_w1,std_w2)
% refer to Shuster M D, S. D. O H. Three-axis attitude determination from vector observations[J]. J Guidance Control Dynam, 1981, 4(1).
var1=std_w1^2;
var2=std_w2^2;
P_phi=var1*eye(3)+1/norm(cross(w1,w2))^2*((var2-var1)*w1*w1'...
    +var1*(w1'*w2)*(w1*w2'+w2*w1'));
RMSE_phi=sqrt(diag(P_phi));