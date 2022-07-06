function clusterPlot(Y, z)

k = length(unique(z));
reLabelTab = table((1:k)', sort(unique(z))','VariableNames',{'cIdx' 'z'});
z = join(table(z','VariableNames',{'z'}), reLabelTab).cIdx;

N = size(Y, 1);
T = size(Y, 2);
YPlot = Y + repmat(linspace(10*N, 10, N)', 1, T);

colTable = table((1:k)', table(hsv(k)),'VariableNames',{'cIdx' 'col'});
clusCol = join(table(z,'VariableNames',{'cIdx'}),colTable);
h = plot(YPlot');
set(h, {'color'}, num2cell(table2array(clusCol.col),2));
set(gca,'ytick',[])

 
 
end