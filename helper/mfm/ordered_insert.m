function actList = ordered_insert(index, actList, t)
j= t;
while((j>0) && (actList(j)>index))
   actList(j+1) = actList(j);
   j = j - 1;
end
actList(j+1) = index;

end