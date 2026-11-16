function [a_candidates,b_candidates,exc,exc_value,epsilon,epsilon2] = ...
    constructMiadCandidates( ...
    i,j,a_i1,a_i2,b_j,vector_angle,sigma_a_i1,sigma_a_i2,sigma_b_j)
% Construct admissible (a,b) candidates from the three primary components.
exc=[0,0,0];
exc_value=[0,0,0];
epsilon=0;
epsilon2=0; %used when i==j
b_j_original=b_j;
[s_ai,a_i1,a_i2,noise_ai,exception,ai1ai2] = ...
    inferMissingComponentCandidatesWithNoise( ...
    a_i1,a_i2,sigma_a_i1,sigma_a_i2);
exc(1)=exception;
exc_value(1)=ai1ai2;
a_candidates=[];
b_candidates=[];
s_a_special=[];
s_b_special=[];
s_change=zeros(2,2);
s_change_count=0;
ctheta=cos(vector_angle);
stheta=sin(vector_angle);
for k=1:1:length(s_ai)
    a=zeros(3,1);
    a(i)=s_ai(k);
    a(cyclicIndex(i,1))=a_i1;
    a(cyclicIndex(i,2))=a_i2;
    aj=a(j);
    aj1=a(cyclicIndex(j,1));
    aj2=a(cyclicIndex(j,2));
    b_j=b_j_original;
    delta=stheta^2-aj^2-b_j^2+2*aj*b_j*ctheta;  
    exc_value(k+1)=delta;
    if delta>0
        [a_candidates,b_candidates]=appendNonzeroDeltaCandidates( ...
            j,aj,aj1,aj2,b_j,delta,ctheta,a_candidates,b_candidates);
    elseif delta==0
            [a_candidates,b_candidates]=appendZeroDeltaCandidates( ...
                j,aj,aj1,aj2,b_j,ctheta,a_candidates,b_candidates);
    else
        switch j
                case i
                    noise_aj=noise_ai;
                    [aj,b_j,is_exist,change]=rationalizeData( ...
                        aj,b_j,vector_angle,noise_aj,sigma_b_j);
                    aj1=aj1*sqrt(1-aj^2)/sqrt(aj1^2+aj2^2);
                    aj2=aj2*sqrt(1-aj^2)/sqrt(aj1^2+aj2^2);
                    if is_exist                        
                        [a_candidates,b_candidates]=appendZeroDeltaCandidates( ...
                            j,aj,aj1,aj2,b_j,ctheta,a_candidates,b_candidates);
                        if k==1
                            epsilon=change;
                            exc(2)=1;
                        else
                            epsilon2=change;
                            exc(3)=1;
                        end
                    else
                        [s_a_special,s_b_special]=appendZeroDeltaCandidates( ...
                            j,aj,aj1,aj2,b_j,ctheta,s_a_special,s_b_special);
                        s_change_count=s_change_count+1;
                        s_change(s_change_count,:)=[k,change];
                    end           
                case cyclicIndex(i,1)
                    noise_aj=sigma_a_i1;
                    [aj,b_j,~,change]=rationalizeData( ...
                        aj,b_j,vector_angle,noise_aj,sigma_b_j);
                    a_i1=aj;
                    [s_ai,a_i1,a_i2]=inferMissingComponentCandidates( ...
                        a_i1,a_i2);
                    aj=a_i1;
                    aj1=a_i2;
                    for n=1:1:length(s_ai)
                        aj2=s_ai(n);
                        [a_candidates,b_candidates]=appendZeroDeltaCandidates( ...
                            j,aj,aj1,aj2,b_j,ctheta,a_candidates,b_candidates);
                    end
                    exc(2)=1;
                    epsilon=change;
                    break;
                case cyclicIndex(i,2)
                    noise_aj=sigma_a_i2;
                    [aj,b_j,~,change]=rationalizeData( ...
                        aj,b_j,vector_angle,noise_aj,sigma_b_j);
                    a_i2=aj;
                    [s_ai,a_i1,a_i2]=inferMissingComponentCandidates( ...
                        a_i1,a_i2);
                    aj=a_i2;
                    aj2=a_i1;
                    for n=1:1:length(s_ai)
                        aj1=s_ai(n);
                        [a_candidates,b_candidates]=appendZeroDeltaCandidates( ...
                            j,aj,aj1,aj2,b_j,ctheta,a_candidates,b_candidates);
                    end
                    exc(2)=1;
                    epsilon=change;
                    break;
        end
    end
end
if isempty(a_candidates)
    if s_change_count == 0
        error('MIAD:NoFeasibleCandidate', ...
            'No feasible attitude candidate could be constructed.');
    end
    [minimum_change,I]=min(s_change(1:s_change_count,2));
    if s_change(I,1)==1
        epsilon=minimum_change;
        exc(2)=1;
    else
        epsilon2=minimum_change;
        exc(3)=1;
    end
    a_candidates=s_a_special(:,I);
    b_candidates=s_b_special(:,I);
end

function [aj,b_j,is_exist,change]=rationalizeData( ...
    aj,b_j,vector_angle,noise_aj,sigma_b_j)
k=noise_aj/sigma_b_j;
ctheta=cos(vector_angle);
diff_aj=-2*aj+2*b_j*ctheta;
diff_bj=-2*b_j+2*aj*ctheta;
if (diff_aj*diff_bj)>0
    a=-1-k^2+2*k*ctheta;
    b=-2*aj*k-2*b_j+2*aj*ctheta+2*b_j*k*ctheta; 
    c=sin(vector_angle)^2-aj^2-b_j^2+2*aj*b_j*ctheta;    
    direction=true;
else
    a=-1-k^2-2*k*ctheta;
    b=2*aj*k-2*b_j+2*aj*ctheta-2*b_j*k*ctheta;
    c=sin(vector_angle)^2-aj^2-b_j^2+2*aj*b_j*ctheta;
    direction=false;
end
sqrt_delta=sqrt(b^2-4*a*c);
x=[(-b+sqrt_delta)/2/a,(-b-sqrt_delta)/2/a];
[~,I]=min(abs(x));
if direction
    aj=aj+k*x(I);
    b_j=b_j+x(I);   
else
    aj=aj-k*x(I);
    b_j=b_j+x(I);
end
if abs(x(I))<4*sigma_b_j    %
    is_exist=true;
else
    is_exist=false;
end
change=x(I);

function [s_ai,a_i1,a_i2,noise_ai,exception,quad_ai1ai2] = ...
    inferMissingComponentCandidatesWithNoise( ...
    a_i1,a_i2,sigma_a_i1,sigma_a_i2)
quad_ai1ai2=a_i1^2+a_i2^2;
exception=0;
if quad_ai1ai2>=1
    exception=1;
    norm_ai1ai2=sqrt(quad_ai1ai2);
    a_i1=a_i1/norm_ai1ai2;
    a_i2=a_i2/norm_ai1ai2;
    s_ai=0;
    noise_max=max([sigma_a_i1,sigma_a_i2]);
    noise_ai=sqrt(sqrt(2))/sqrt(noise_max)*noise_max;   % Cap propagated noise using the feasibility-bound threshold.
else
    s_ai=sqrt(1-quad_ai1ai2)*[1,-1];
    noise_ai=sqrt((a_i1/s_ai(1))^2*sigma_a_i1^2+(a_i2/s_ai(1))^2*sigma_a_i2^2);
    noise_max=max([sigma_a_i1,sigma_a_i2]);
    noise_ai=min(noise_ai,sqrt(2*sqrt(2))/sqrt(noise_max)*noise_max);
end

function [s_ai,a_i1,a_i2]=inferMissingComponentCandidates(a_i1,a_i2)
quad_ai1ai2=a_i1^2+a_i2^2;
if quad_ai1ai2>=1
    norm_ai1ai2=norm([a_i1,a_i2]);
    a_i1=a_i1/norm_ai1ai2;
    a_i2=a_i2/norm_ai1ai2;
    s_ai=0;
else
    s_ai=sqrt(1-quad_ai1ai2)*[1,-1];
end

function [a_candidates,b_candidates]=appendZeroDeltaCandidates( ...
    j,aj,aj1,aj2,b_j,ctheta,a_candidates,b_candidates)
quad_aj1aj2=aj1^2+aj2^2;
bj1=aj1*(ctheta-aj*b_j)/quad_aj1aj2;
bj2=aj2*(ctheta-aj*b_j)/quad_aj1aj2;
a=zeros(3,1);
a(j)=aj;
a(cyclicIndex(j,1))=aj1;
a(cyclicIndex(j,2))=aj2;
b=zeros(3,1);
b(j)=b_j;
b(cyclicIndex(j,1))=bj1;
b(cyclicIndex(j,2))=bj2;
a_candidates=[a_candidates,a];
b_candidates=[b_candidates,b];

function [a_candidates,b_candidates]=appendNonzeroDeltaCandidates( ...
    j,aj,aj1,aj2,b_j,delta,ctheta,a_candidates,b_candidates)
quad_aj1aj2=aj1^2+aj2^2;
sqrt_delta=sqrt(delta);
s_bj1=(aj1*(ctheta-aj*b_j)+[1,-1]*aj2*sqrt_delta)/quad_aj1aj2;
s_bj2=(aj2*(ctheta-aj*b_j)+[-1,1]*aj1*sqrt_delta)/quad_aj1aj2;
a=zeros(3,1);
a(j)=aj;
a(cyclicIndex(j,1))=aj1;
a(cyclicIndex(j,2))=aj2;
b=zeros(3,1);
b(j)=b_j;
b(cyclicIndex(j,1))=s_bj1(1);
b(cyclicIndex(j,2))=s_bj2(1);
a_candidates=[a_candidates,a];
b_candidates=[b_candidates,b];
b=zeros(3,1);
b(j)=b_j;
b(cyclicIndex(j,1))=s_bj1(2);
b(cyclicIndex(j,2))=s_bj2(2);
a_candidates=[a_candidates,a];
b_candidates=[b_candidates,b];