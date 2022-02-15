function MultiphaseDisp(varargin)
%% image setting
TitleSettings = {[],'FontSize',20,'Fontname','arial'};%,'fontweight','bold'};
sigma = [];
ContactComp = [];
color = [];
FaceAlpha = 0.7;
th = 0.35;
PhiStack = varargin{1};
CellIdx = varargin{2};
[CellIdx,I] = sort(CellIdx);
PhiStack = PhiStack(I);
dr = varargin{3};
CellParameters = varargin{4};
TripleView = 0;
EnableLegend = 1;
EnableBox = 0;
N = size(PhiStack{1});
L = [zeros(1,3),N([2,1,3])*dr];
AxisLimit = [];
ViewVec = [0,0,-1];
for StrIdx = find(cellfun(@ischar,varargin))
	switch varargin{StrIdx}
		case 'title'
			TitleSettings{1} = varargin{StrIdx+1};
		case 'threshold'
			th = varargin{StrIdx+1};
		case 'sigma'
			sigma = varargin{StrIdx+1};
      sigma = sigma(I,I);
		case 'FaceAlpha'
			FaceAlpha = varargin{StrIdx+1};
		case 'TripleView'
			TripleView = true;
    case 'ContactComp'
      ContactComp = varargin{StrIdx+1};
    case 'LegendOff'
      EnableLegend=0;
    case 'color'
      color = varargin{StrIdx+1}(I,:);
    case 'AxisLimit'
      AxisLimit = varargin{StrIdx+1};
    case 'DomainRange'
      L = varargin{StrIdx+1};
    case 'view'
      ViewVec = varargin{StrIdx+1};
    case 'AddBox'
      EnableBox = 1;
	end
end

%% config cell name, color
if any(contains(fieldnames(CellParameters),'FaceAlpha'))
  FaceAlphaList = CellParameters.FaceAlpha(CellIdx);
else
  FaceAlphaList = FaceAlpha*ones(length(CellIdx),1);
end
if isempty(color)
  color = CellParameters.Color(CellIdx,:);
end
if isempty(AxisLimit)
  AxisLimit = L;
end
if EnableLegend
  CellName = CellParameters.Name(CellIdx);
  LegendSettings = {CellName,'FontSize',17.5,'FontName','arial','box','off'};
  FontSize = 20;
else
  FontSize = 20;%25.5;
end

%% iso-surface
[X,Y,Z] = meshgrid(L(1):dr:L(4)-dr,L(2):dr:L(5)-dr,L(3):dr:L(6)-dr);
fv = cellfun(@(x)isosurface(X,Y,Z,x,th),PhiStack,'UniformOutput',0);

%% imaging
if TripleView
	y0 = 0.2;
	x0 = 0.1;
	xw = [15,30,30]*0.75/45;
	yw = [20,20,15]*0.75/45;
	ViewDir = [
    -1,0,0;% DV-LR
    0,0,-1;% AP-DV
    0,-1,0];% AP-LR
  CamRollAngle = [90,0,0];
	ObjStack = cell(3,1);
  for ViewIdx = 1:3
    ObjStack{ViewIdx} = axes('Position',[x0,y0,xw(ViewIdx),yw(ViewIdx)]);
    IsoSurfacePlot(fv,color,FaceAlphaList,20,ViewDir(ViewIdx,:),L,X,Y,Z,PhiStack);
    camroll(CamRollAngle(ViewIdx))
    if ViewIdx<3
      if xw(ViewIdx) ~= xw(ViewIdx+1)
        x0 = x0+0.1+xw(ViewIdx);
      end
      if yw(ViewIdx) ~= yw(ViewIdx+1)
        y0 = y0+0.1+yw(ViewIdx);
      end
    end
  end
  set(ObjStack{1},'YAxisLocation','right')
	if ~isempty(sigma)
		% legend, adjacent matrix and sigma
    ax = {axes('Position',[0.1,y0,0.25,0.25]),axes('Position',[0.01,y0,0.08*0.25,0.25])};
  else
    % only legend
		h = findobj(gca);
		legend(h(end:-1:3),LegendSettings{:},'position',[0,0.76,0.425,0]);
	end
	sgtitle(TitleSettings{:});
else
	width = 0.32;
	gap = 0.05;
	y0 = 0.5-width/2;
	if isempty(sigma)
    title(TitleSettings{1:5});
  else
		axes('Position',[width+gap*3.2,y0,1.5*width,width]);
  end
  IsoSurfacePlot(fv,color,FaceAlphaList,FontSize,ViewVec,AxisLimit,X,Y,Z,PhiStack);
  if EnableBox
    AddBox(AxisLimit)
  end
  if EnableLegend
    if isempty(sigma)
      h = findobj(gca);
      legend(h(end:-1:3),LegendSettings{:},'Location','eastoutside','NumColumns',1+(length(CellIdx)>26));
    else
      % adjacent matrix and sigma
      sgtitle(TitleSettings{:});
      ax = {axes('Position',[1.7*gap,y0,width,width]),axes('Position',[0.01,y0,0.08*width,width])};
    end
  end
end
if ~isempty(sigma)&&EnableLegend
  AdjMatSigmaLegend(IsContact(PhiStack,0.44),ContactComp,sigma,CellName,ax,color)
end
end

%% subfunctions
function PatchStack = IsoSurfacePlot(fv,color,FaceAlphaList,FontSize,ViewVec,AxisLimit,X,Y,Z,PhiStack)
CellNum = length(fv);
PatchStack = cell(CellNum,1);
for n = 1:CellNum
  PatchStack{n} = patch(fv{n},'FaceColor',color(n,:),'EdgeColor','none','FaceAlpha',FaceAlphaList(n));
  isonormals(X,Y,Z,PhiStack{n},PatchStack{n})
end
axis('equal',AxisLimit([1,4,2,5,3,6]))
grid off
ylabel({'\itz\rm / D-V Axis (\mum)'});
xlabel({'\itx\rm / A-P Axis (\mum)'});
zlabel({'\ity\rm / L-R Axis (\mum)'});
set(gca,'FontSize',FontSize,'Fontname','arial','ydir','reverse','tickdir','out','LineWidth',2);
view(ViewVec)
camlight(-20,40)
lighting gouraud
end

function AddBox(L)
hold on
[x,y,z] = meshgrid(L([1,4]),L([2,5]),L([3,6]));
v = [x(:),y(:),z(:)];
plotted = zeros(8,1);
for n = 1:8
  idx = find(sum(v(n,:)==v,2)==2);
  for m = find(plotted(idx)<3)'
    plot3(v([n,idx(m)],1),v([n,idx(m)],2),v([n,idx(m)],3),'k');
  end
  plotted(n) = 3;
  plotted(idx) = plotted(idx)+1;
end
end
% function AddBackground(L)
% [x,y,z] = meshgrid(L([1,4]),L([2,5]),L([3,6]));
% v = [x(:),y(:),z(:)];
% f = [1,3,7,5;1,2,4,3;3,4,8,7];
% patch('Faces',f,'Vertices',v,'FaceColor','w','EdgeColor','none','FaceAlpha',0.5)
% end