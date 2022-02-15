function HighlightSquare(data,color,ax)
if nargin<=2
  ax = gca;
end
CellNum = size(data,1);
f = 1:4;
[i,j] = find(data);
for n = 1:length(i)
  v = [i(n)-[0;1;1;0],j(n)-[0;0;1;1]];
  patch(ax,'Faces',f,'Vertices',v,'FaceColor','none','EdgeColor',color,'LineWidth',2*18/CellNum)
end