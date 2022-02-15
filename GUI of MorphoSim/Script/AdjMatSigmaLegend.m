function AdjMatSigmaLegend(DispContact,CompContact,sigma,CellName,ax,CellColor)
if nargin<=4
  ax{1} = gca;
end
CellNum = length(CellName);
% sigma
if sigma == 0
  sigma = ones(CellNum);
  SigmaValue = 1;
  c = 1;
else
  SigmaValue = unique(sigma(:));
  if length(SigmaValue)>2
    c=[1,0.8*2.^-(0:length(SigmaValue)-2)]';
  else
    c = [1;0.4];
  end
end
SquarePlot(sigma,CellName,[SigmaValue,c*ones(1,3)],ax{1})
% adjacency matrix 
MarkerSetting = {1, 'k.', 180/CellNum;% conserve
                 2, 'ko', 180/CellNum/4};%non-conserve
SquarePlotMarker(DispContact,MarkerSetting,ax{1})
% highlight broken conserved contact
if ~isempty(CompContact)
  if all(size(DispContact) == size(CompContact))
    if length(unique(CompContact))>2
      SimContact = DispContact;
      ExpContact = CompContact;
    else
      ExpContact = DispContact;
      SimContact = CompContact;
    end
    HighlightSquare(((SimContact-(ExpContact==1))>0&ExpContact~=2)|(SimContact-(ExpContact==1))<0,[1,0,0],ax{1});
  else
    warning('Contact map size is not consistant')
  end
end
% legend
if nargin > 5
  f = 1:4;
  w = 3*6/max([3,CellNum]);
  for n = 1:CellNum
    vert = [[0;1;1;0],[0;0;1;1]-n];
    patch(ax{2},'Faces',f,'Vertices',vert,'FaceColor',CellColor(n,:),'EdgeColor',[1,1,1],'LineWidth',w)
  end
  axis(ax{2},[0,1,-CellNum,0])
  set(ax{2},'xticklabels',[],'yticklabels',[])
end