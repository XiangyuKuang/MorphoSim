function SquarePlot(data,CellName,color,ax)
% color = [1,  0,0,1;
%          0,  0.5,0.5,0.5;
%         -1,  1,0,0;
%         -2,  1,1,1]
if nargin<=3
  ax = gca;
end
CellNum = length(CellName);
f = [1:4;5:8];
for n = 1:CellNum
  for m = n:CellNum 
    v = [[n-[0;1;1;0],m-[0;0;1;1]];[m-[0;1;1;0],n-[0;0;1;1]]];
    patch(ax,'Faces',f,'Vertices',v,'FaceColor',color(color(:,1) == data(n,m),2:4),'EdgeColor',ones(1,3),'LineWidth',18/CellNum)
  end
end
dx = -1.2/CellNum;
if CellNum>8
  fontsize = 144/CellNum;
else
  fontsize = 15;
end
axis(ax,'equal',[0,CellNum,0,CellNum])
set(ax,'FontSize',fontsize,'Fontname','arial','ydir','reverse',...
  'xtick',0.5+dx:1:0.5+CellNum+dx,'xticklabel',CellName,'xticklabelrotation',60,...
  'ytick',0.5:1:CellNum+0.5,'yticklabel',CellName)