function [ album,swap ] = updatealbum2( album,swap )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
    booster=draw(8,640)';
    album(booster)=album(booster)+1;
    swap=swap + (album>1);
    album=double(album&album);
end