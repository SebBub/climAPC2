function tmap=sigtrendmap(M,time)

% TMAP=SIGTRENDMAP(M);
%
% M: Input monthly data matrix <time x lat x lon>
% time: time vector
%
% tmap: Map with 95% significant trends <lat x lon> in units of data units
% over time units (e.g. centigrade per day if the time vector is in days)

X=ones(size(M,1),2);
%X(:,2)=[1:size(M,1)]';
X(:,2)=time;

tmap=zeros(size(M,2),size(M,3));

id='stats:regress:RankDefDesignMat';
warning('off',id)

for y=1:size(M,2)
    for x=1:size(M,3)
        yM=M(:,y,x);
        [b,bint]=regress(yM,X);
        tmap(y,x)=b(2);
        pm=(bint(2,2)-bint(2,1))/2;
        if ( abs(tmap(y,x)) <= abs(pm) )
            tmap(y,x)=NaN;
        end
        tmap(y,x)=tmap(y,x);
    end
end

