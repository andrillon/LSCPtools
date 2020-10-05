function v2=minmax(v,dim)

if nargin<2
    dim=2;
end

if min(size(v))==1
    v2=(v-nanmin(v))/nanmax(v-nanmin(v));
else
    if dim==1
        v2=(v-repmat(nanmin(v,[],dim),size(v,1),1))./repmat(nanmax(v-repmat(nanmin(v,[],dim),size(v,1),1),[],dim),size(v,1),1);
    elseif dim==2
        v2=(v-repmat(nanmin(v,[],dim),1,size(v,2)))./repmat(nanmax(v-repmat(nanmin(v,[],dim),1,size(v,2)),[],dim),1,size(v,2));
    end
end