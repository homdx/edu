function [ album,swap ] = createAlbum2( q )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
n=640;
album=zeros(1,n);
swap=zeros(1,n);
for qq=1:q
    booster=draw(8,640)';
    album(booster)=album(booster)+1;
    swap=swap + (album>1);
    album=double(album&album);
end
end

