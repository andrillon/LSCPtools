function [closestvalue,index]=findclosest(vec,value)

if length(value)==1
    vecdiff=vec-value;
    index0=find(abs(vecdiff)==min(abs(vecdiff)));
    if isempty(index0)
        closestvalue=[];
        index=[];
    else
        index=index0(1);
        closestvalue=vec(index0(1));
    end
    
else
    for k=1:length(value)
        vecdiff=vec-value(k);
        index0=find(abs(vecdiff)==min(abs(vecdiff)));
        if isempty(index0)
            closestvalue(k)=NaN;
            index(k)=NaN;
        else
            index(k)=index0(1);
            closestvalue(k)=vec(index0(1));
        end
    end
end