function RMSE_phi=estimateQuestRmse(w1,w2,std_w1,std_w2)
% refer to Shuster M D, S. D. O H. Three-axis attitude determination from vector observations[J]. J Guidance Control Dynam, 1981, 4(1).

%%%%%%% theoretical result in the Shuster's paper
var1=std_w1^2;
var2=std_w2^2;
var_tot=1/(1/var1+1/var2);
P_phi=var_tot*eye(3)+(norm(cross(w1,w2)))^(-2)*((var2-var_tot)*w1*w1'...
    +(var1-var_tot)*w2*w2'+var_tot*(w1'*w2)*(w1*w2'+w2*w1'));
RMSE_phi=sqrt(diag(P_phi));