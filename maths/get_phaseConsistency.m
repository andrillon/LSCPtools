function PC=get_phaseConsistency(phase)

if ndims(phase)<3
    if min(size(phase))==1
        PC=sqrt((1/length(phase)*sum(sin(phase))).^2+(1/length(phase)*sum(cos(phase))).^2);
    else
        for n=1:size(phase,1)
            PC(n,:)=sqrt((1/length(phase(n,:))*sum(sin(phase(n,:)))).^2+(1/length(phase(n,:))*sum(cos(phase(n,:)))).^2);
        end
    end
else
   for i=1:size(phase,1)
       for j=1:size(phase,2)
           temp=squeeze(phase(i,j,:));
           PC(i,j)=sqrt((1/length(temp)*sum(sin(temp))).^2+(1/length(temp)*sum(cos(temp))).^2);
       end
   end    
end