function v=extractDiscriminatorComponents(a,b,m,n)
[row,col]=size(m);
if row==1
    m=m';
    n=n';
end
[row,col]=size(m);
v=zeros(3,col);
v(1,:)=m(a,:);
v(2,:)=n(cyclicIndex(b,1),:);
v(3,:)=n(cyclicIndex(b,2),:);