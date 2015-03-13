clear all;
clc;
stickeramount=0;
n=input('Number of collectors (Warning: high numbers may result in high run time): ');
rep=input('Number of iterations(Warning in creases runtime linearly): ');
startpacknumber=input('Number of packages bought before starting to swap: ');
AlbCell=cell(n,2);
albsum=zeros(n,1);

tic
for q=1:rep;
    x=1:n;
    % [alb1 swap1]=updatealbum([],[]);
    % [alb2 swap2]=updatealbum([],[]);
    for t=1:n
        [AlbCell{t,1} AlbCell{t,2}]=createAlbum2(startpacknumber);
        albsum(t)=sum(AlbCell{t,1});
    end
    while (min(albsum)<640)
        %swapping algorithm here
        if max(albsum)==640
            x=find(albsum~=640)';
        end
        count=1;
        for pp=x(1:length(x)-1)
            count=count+1;
            y=x(count:(length(x)));
            for qq=y
                [AlbCell{pp,1} AlbCell{pp,2} AlbCell{qq,1} AlbCell{qq,2}]=swap2(AlbCell{pp,1},AlbCell{pp,2},AlbCell{qq,1},AlbCell{qq,2});
            end
        end
        for ll=x
            for nn=1:5
            [AlbCell{ll,1} AlbCell{ll,2}]=updatealbum2(AlbCell{ll,1},AlbCell{ll,2});
            end;
            albsum(ll)=sum(AlbCell{ll,1});
        end
    end
    stickeramount=sum(albsum)+stickeramount;
    for rr=1:n
        stickeramount=stickeramount+sum(AlbCell{rr,2});
    end
end
stickeramount=stickeramount/rep
toc
