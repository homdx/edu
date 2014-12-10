function [ sticker ] = draw( packsize,numofstickers )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
n=packsize;
num=numofstickers;
if n<=num
cards=1+floor(rand(n,1)*num);
i=1;
j=1;
p=0;
    while (p==0)
        for i=1:n-1
            for j=i+1:n
                if cards(i)==cards(j)
                    cards(i)=1+floor(rand(1)*num);
                    p=-1;
                end
            end
        end
        if p==-1
            p=0;
        else p=1;
        end
    end
sticker=cards;
else
    sticker=[];
end
end

