function neighber_index=find_neighber(ind_min,tri)
% find the neighber
neighber_index=[];
for ii=1:size(tri,1)
  if ismember(ind_min,tri(ii,:))
    neighber_index=union(neighber_index,tri(ii,:));
  end
end
neighber_index(neighber_index==ind_min)=[];
end