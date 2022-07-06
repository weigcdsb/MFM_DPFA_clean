function actList = ordered_remove(index, actList, t)
for j = 1:t
   if(actList(j) >= index)
      actList(j) = actList(j+1); 
   end
end
end
