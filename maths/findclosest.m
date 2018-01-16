function [closestvalue,index]=findclosest(vec,value)

vecdiff=vec-value;
index=find(abs(vecdiff)==min(abs(vecdiff)));
if isempty(index)
    closestvalue=[];
else
    index=index(1);
    closestvalue=vec(index(1));
end