function [ alb1,swap1,alb2,swap2 ] = swap2( alb1,swap1,alb2,swap2 )
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here
A1=~alb1&swap2;
A2=~alb2&swap1;
ind1=[];
ind2=[];
if 0<sum(A1) && sum(A1)<=sum(A2)
    ind1=find(A1);
    ind2=find(A2,sum(A1),'first');
end
if 0<sum(A2)&& sum(A2)<sum(A2)
    ind2=find(A2);
    ind1=find(A1,sum(A1),'first');
end
alb1(ind1)=alb1(ind1)+1;
alb2(ind2)=alb2(ind2)+1;
swap1(ind2)=swap1(ind2)-1;
swap2(ind1)=swap2(ind1)-1;

end