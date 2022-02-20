function SquarePlotMarker(data,MarkerSetting,ax)
if nargin<=2
  ax = gca;
end
hold(ax,'on');
for n = 1:size(MarkerSetting,1)
  [y,x] = find(data==MarkerSetting{n,1});
  plot(ax,x-0.5,y-0.5,MarkerSetting{n,2},'markersize',MarkerSetting{n,3})
end