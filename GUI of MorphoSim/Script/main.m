function SimReport = main(SimName,L,N,EggshellPara,h,SimParameters,CellParameters,StageEnd,TimeSpan,varargin)
% BDF2 scheme
% SimReport = {N,h,StageEndTime,SimEndInfo,escapeT,trajectory,velocity,CellIdx};
%   StageEndTime format: 2x(N-1), N = StageEndNum
%             |StageEnd(2)|,...,|StageEnd(N)|
%            T|           |     |           |
%     TimeCost|           |     |           |
% difference between non-gui version and this code: 
% 1) add EnableEggshell at line 35 and 72
% 2) print logs in Commond Window instead of a log file 
% 3) disable saving StageEnd.mat at line 73
% 4) disable cell division ar line 118
% 5) disable SimParaModify at line 126
% 6) convert gamma to vector at line 129
% 7) disable saving initial condition at line 156
% 8) disable saving trajectory.mat and velocity.mat at line 306
% 9) auto delete CellParameter.mat at line 307
% 
% Xiangyu Kuang
%% input config
dr = L./N;
dv = prod(dr);
[xx,yy,zz] = meshgrid(0:dr(1):L(1)-dr(1),0:dr(2):L(2)-dr(2),0:dr(3):L(3)-dr(3));% periodic, phi(0)=phi(L)
GPU = false;
SavePhi = false;
LoadDir = SimName;
s = [0,0];
LoadSwitch = TimeSpan(1) > 0;
EnforceStageIdx = [];
for StrIdx = find(cellfun(@ischar,varargin))
  switch varargin{StrIdx}
    case 'SavePhi'
      SavePhi = true;
      SaveStep = varargin{StrIdx+1};
    case 'Uncompressed'
      EggshellPara(4) = -1;
    case 'GPU'
      GPU = true;
    case 'LoadDir'
      LoadSwitch = true;
      LoadDir = varargin{StrIdx+1};
    case 'StabilizationPara'
      s = varargin{StrIdx+1};
    case 'EnforceStage'
      EnforceStageIdx = varargin{StrIdx+1};
    case 'EnableEggshell'
      EnableEggshell = varargin{StrIdx+1};
  end
end
if exist(SimName,'dir') ~= 7% make dir named by SimName
  mkdir(SimName)
end
if GPU && SavePhi% init parpool for asynchronously saveing(ParaSave)
  p = gcp;
  p.IdleTimeout = inf;
  f = parfeval(p,@rand,1);
end
%% init log
% [fileID, msg] = fopen(sprintf('%s/%s_log',SimName,SimName),'a');
%% load or generate initial state
if LoadSwitch% continue running
  % load initial state
  [PhiStack,CellIdx,StageBreak,StageEndIdx] = LoadInitState(TimeSpan(1),LoadDir,StageEnd);
else % new simulation, never used in gui
  [PhiStack,CellIdx,StageBreak,StageEndIdx] = GenerateInitState(CellParameters,EggshellPara,xx,yy,zz,L,dv,StageEnd);
end
if ~isempty(EnforceStageIdx)
  StageEndIdx = EnforceStageIdx:length(StageEnd)-1;
  if isempty(StageEndIdx)
    error('Enforce StageIdx larger than length(StageEnd)-1')
  end
end
% save StageEnd % removed
CellNum = length(CellIdx);
% load or create AnalVar
trajectory = InitAnalVar(SimName,'trajectory',LoadSwitch);
velocity = InitAnalVar(SimName,'velocity',LoadSwitch);
% generate eggshell
if EnableEggshell
  PhiShell = (1-tanh(ShapeInit(EggshellPara,L,xx,yy,zz)/eps))/2;
  EggshellChecking(PhiStack,PhiShell,0.34,N)
else
  PhiShell = 0;
end

%% init var
% wavenum
[kx,ky,kz]=meshgrid(fftshift(-floor(N(1)/2):ceil(N(1)/2)-1),fftshift(-floor(N(2)/2):ceil(N(2)/2)-1),fftshift(-floor(N(3)/2):ceil(N(3)/2)-1));
kx = 2*pi/L(1)*1i*kx;
ky = 2*pi/L(2)*1i*ky;
kz = 2*pi/L(3)*1i*kz;
% transfer var to gpu if GPU is true
if GPU
  PhiStack = cellfun(@gpuArray,PhiStack,'UniformOutput',0);
  PhiShell = gpuArray(PhiShell);
  kx = gpuArray(kx);
  ky = gpuArray(ky);
  kz = gpuArray(kz);
  xx = gpuArray(xx);
  yy = gpuArray(yy);
  zz = gpuArray(zz);
end
ksq = kx.^2+ky.^2+kz.^2;
fPhiStack = cellfun(@fftn,PhiStack,'UniformOutput',0);
t = TimeSpan(1);% time in simulation
TimeCost = 0;% runing time
SimBreak = false;% terminate simulation if true, updated by IsEnd
AnalStep = floor(10/h);% do analysis every AnalStep
StageEndTime = zeros(2,StageEndIdx(end));% record t when StageEnd

%% simulation
fprintf('Simulation %s in progress...\n%0.1f Init\n',SimName,t);
tic
for StageIdx = StageEndIdx
  if SimBreak% terminate simulation if true
    continue
  end
  % cell division % removed
  % init velocity&trajectory
  r = cell(CellNum+1,1);
  v = cell(CellNum+1,1);
  for n = 1:CellNum
    r{n} = gather([sum(PhiStack{n}.*xx,'all'),sum(PhiStack{n}.*yy,'all'),sum(PhiStack{n}.*zz,'all')]/fPhiStack{n}(1));
  end
  r{end} = t;
  % update parameters
  tau = SimParameters.tau;
  M = SimParameters.M;
  gamma = SimParameters.gamma*ones(CellNum,1);
  c = SimParameters.c;
  g = SimParameters.g;
  g_shell = SimParameters.g_shell;
  sigma = SimParameters.sigma;
  Volume = CellParameters.Volume(CellIdx);
  MR = M./Volume;
  % init sim var
  [gammaValue,~,gammaIdx] = unique(gamma);
  sigma1 = zeros(1,CellNum);
  if length(unique(sigma(:))) == 1
    sigma1 = sigma(1)*ones(1,CellNum);
    sigmaIdx = zeros(CellNum);
  else
    for n = 1:CellNum
      sigma1(n) = mode(sigma(sigma(:,n)~=0,n));
    end
    sigmaIdx = (sigma~=sigma1)-eye(CellNum);
  end
  sigmaSub = sigma1-sigma;
  EggShellRep = g_shell*PhiShell.^2;
  % init stage end
  StageBreak = false;
  NextStageEnd = StageEnd(StageIdx+1).StageEndDef;
  [DefSwitch, DefIdx]= ismember({'Duration';'SS';'QSS'},NextStageEnd(:,1));
  t0 = t;
  countQSS = 0;% count QSS for StageEndDef
  % save IC % removed
  % computation
  PhiXStack = cell(CellNum,1);
  PhiYStack = cell(CellNum,1);
  PhiZStack = cell(CellNum,1);
  PhiSq = cell(CellNum,1);
  fFOld = cell(CellNum,1);
  fPhiStackOld = fPhiStack;
  A = arrayfun(@(x)(tau/h+s(1)-x*ksq).^-1,gammaValue,'UniformOutput',0);
  B = tau/h+s(1);
  C = 0.5*tau/h+s(1);
  D = 2*gamma*c;
  % compute phi^1_i with first order scheme
  tic
  IntervalLength = AnalStep-1;% IntervalLength will be updated by IsEnd
  PhiXStack{1} = real(ifftn(fPhiStack{1}.*kx));
  PhiYStack{1} = real(ifftn(fPhiStack{1}.*ky));
  PhiZStack{1} = real(ifftn(fPhiStack{1}.*kz));
  PhiSq{1} = PhiStack{1}.^2;
  Rep = PhiSq{1};
  PhiXSum = PhiXStack{1};
  PhiYSum = PhiYStack{1};
  PhiZSum = PhiZStack{1};
  for n = 2:CellNum
    PhiXStack{n} = real(ifftn(fPhiStack{n}.*kx));
    PhiYStack{n} = real(ifftn(fPhiStack{n}.*ky));
    PhiZStack{n} = real(ifftn(fPhiStack{n}.*kz));
    PhiSq{n} = PhiStack{n}.^2;
    Rep = Rep + PhiSq{n};
    PhiXSum = PhiXSum + PhiXStack{n};
    PhiYSum = PhiYSum + PhiYStack{n};
    PhiZSum = PhiZSum + PhiZStack{n};
  end
  Rep = EggShellRep + g*Rep;
  MRV = MR.*(Volume-cellfun(@(x)real(x(1)),fPhiStack)*dv);
  for n = 1:CellNum
    AtrX = sigma1(n)*(PhiXSum-PhiXStack{n});
    AtrY = sigma1(n)*(PhiYSum-PhiYStack{n});
    AtrZ = sigma1(n)*(PhiZSum-PhiZStack{n});
    for m = find(sigmaIdx(:,n))'
      AtrX = AtrX - sigmaSub(m,n)*PhiXStack{m};
      AtrY = AtrY - sigmaSub(m,n)*PhiYStack{m};
      AtrZ = AtrZ - sigmaSub(m,n)*PhiZStack{m};
    end
    fF = fftn( MRV(n)*sqrt(PhiXStack{n}.^2+PhiYStack{n}.^2+PhiZStack{n}.^2)...%volume constraint
      -(Rep-g*PhiSq{n}).*PhiStack{n}...%repulsion
      -AtrX.*PhiXStack{n}-AtrY.*PhiYStack{n}-AtrZ.*PhiZStack{n}...%attraction
      -D(n)*(PhiSq{n}.*(2*PhiStack{n}-3)+PhiStack{n}));%double well
    fPhiStack{n} = A{gammaIdx(n)}.*(B.*fPhiStack{n}+fF);
    PhiStack{n} = real(ifftn(fPhiStack{n}));
    fFOld{n} = fF;
  end
  TimeCost = TimeCost + toc;
  t = t+h;
  % compute phi^n_i with BDF2 scheme for n>1
  A = arrayfun(@(x)(1.5*tau/h+s(1)-x*ksq).^-1,gammaValue,'UniformOutput',0);
  B = 2*tau/h+2*s(1);
  while ~(StageBreak||SimBreak)
    tic
    for iter = 1:IntervalLength
      PhiXStack{1} = real(ifftn(fPhiStack{1}.*kx));
      PhiYStack{1} = real(ifftn(fPhiStack{1}.*ky));
      PhiZStack{1} = real(ifftn(fPhiStack{1}.*kz));
      PhiSq{1} = PhiStack{1}.^2;
      Rep = PhiSq{1};
      PhiXSum = PhiXStack{1};
      PhiYSum = PhiYStack{1};
      PhiZSum = PhiZStack{1};
      for n = 2:CellNum
        PhiXStack{n} = real(ifftn(fPhiStack{n}.*kx));
        PhiYStack{n} = real(ifftn(fPhiStack{n}.*ky));
        PhiZStack{n} = real(ifftn(fPhiStack{n}.*kz));
        PhiSq{n} = PhiStack{n}.^2;
        Rep = Rep + PhiSq{n};
        PhiXSum = PhiXSum + PhiXStack{n};
        PhiYSum = PhiYSum + PhiYStack{n};
        PhiZSum = PhiZSum + PhiZStack{n};
      end
      Rep = EggShellRep + g*Rep;
      MRV = MR.*(Volume-cellfun(@(x)x(1),fPhiStack)*dv);
      for n = 1:CellNum
        AtrX = sigma1(n)*(PhiXSum-PhiXStack{n});
        AtrY = sigma1(n)*(PhiYSum-PhiYStack{n});
        AtrZ = sigma1(n)*(PhiZSum-PhiZStack{n});
        for m = find(sigmaIdx(:,n))'
          AtrX = AtrX - sigmaSub(m,n)*PhiXStack{m};
          AtrY = AtrY - sigmaSub(m,n)*PhiYStack{m};
          AtrZ = AtrZ - sigmaSub(m,n)*PhiZStack{m};
        end
        fF = fftn( MRV(n)*sqrt(PhiXStack{n}.^2+PhiYStack{n}.^2+PhiZStack{n}.^2)...%volume constraint
          -(Rep-g*PhiSq{n}).*PhiStack{n}...%repulsion
          -AtrX.*PhiXStack{n}-AtrY.*PhiYStack{n}-AtrZ.*PhiZStack{n}...%attraction
          -D(n)*(PhiSq{n}.*(2*PhiStack{n}-3)+PhiStack{n}));%double well
        fPhi_temp = fPhiStack{n};
        fPhiStack{n} = A{gammaIdx(n)}.*(B.*fPhi_temp - C.*fPhiStackOld{n} + 2*fF-fFOld{n});
        PhiStack{n} = real(ifftn(fPhiStack{n}));
        fPhiStackOld{n} = fPhi_temp;
        fFOld{n} = fF;
      end
    end
    TimeCost = TimeCost + toc;
    t = t+IntervalLength*h;
    % analysis
    for n = 1:CellNum
      rc = gather([sum(PhiStack{n}.*xx,'all'),sum(PhiStack{n}.*yy,'all'),sum(PhiStack{n}.*zz,'all')]/real(fPhiStack{n}(1)));
      r{n} = cat(1,r{n},rc);
      v{n} = cat(1,v{n},(r{n}(end,:)-r{n}(end-1,:))/h/IntervalLength);
    end
    r{end} = cat(1,r{end},t);
    v{end} = cat(1,v{end},t);
    % StageEndCheck
    [StageBreak, SimBreak, IntervalLength, StageEndInfo, SimEndInfo, countQSS] = IsEnd(NextStageEnd, DefSwitch, DefIdx, PhiStack, v(1:end-1), countQSS, t, t0, TimeSpan(2), AnalStep, h);
    if DefSwitch(3)
      if strcmp('QSS',StageEndInfo)% return to last analysis step(t-IntervalLength*h) if QSS
        PhiStack = PhiOld;
        fPhiStack = fPhiOld;
        t = tOld;
      end
      PhiOld = PhiStack;
      fPhiOld = fPhiStack;
      tOld = t;
    end
    % save
    if SavePhi && (mod(length(v{end}),SaveStep) == 0||StageBreak||SimBreak)
      if StageBreak
        StageIdxSave = StageIdx+1;
      else
        StageIdxSave = [StageIdx;StageIdx+1];
      end
      if GPU
        % asynchronously saveing
        fetchNext(f);
        f = parfeval(p,@ParaSave,0,'%s/%s_%s.mat',{SimName,SimName,dot2d(t)},cellfun(@gather,PhiStack,'UniformOutput',0),CellIdx,SimParameters,StageIdxSave);
      else
        ParaSave('%s/%s_%s.mat',{SimName,SimName,dot2d(t)},PhiStack,CellIdx,SimParameters,StageIdxSave);
      end
    end
  end
  % update analysis result
  StageEndTime(1,StageIdx) = t;
  StageEndTime(2,StageIdx) = TimeCost;
  trajectory = cat(2,trajectory,{r;CellIdx});
  velocity = cat(2,velocity,{v;CellIdx});
end

%% close
SimReport = {N,h,StageEndTime,SimEndInfo,TimeCost,trajectory,velocity,CellIdx};
if ~StageBreak
  ParaSave('%1$s/%1$s_%2$s.mat',{SimName,dot2d(t)},cellfun(@gather,PhiStack,'UniformOutput',0),CellIdx,SimParameters,StageIdx)
end
% save trajectory and velocity % removed
% delete(sprintf('%1$s/%1$s_CellParameters.mat',SimName));
fprintf('Simulation elapsed time is %0.3fs.\n',TimeCost);
if SavePhi && GPU
  fetchNext(f);
end
end

%% subfunctions
function [PhiStack,CellIdx,StageBreak,StageEndIdx] = GenerateInitState(CellParameters,EggshellPara,xx,yy,zz,L,dv,StageEnd)
StageBreak = true;
CellList = {'P0'};
CellIdx = find(strcmp(CellParameters.Name,CellList));
Volume = CellParameters.Volume(CellIdx);
Rcell = fzero(@(r)Volume-dv*sum((1+tanh(ShapeInit(r*EggshellPara,L,xx,yy,zz)))/2,'all'),1);
PhiStack = {(1+tanh(ShapeInit(Rcell*EggshellPara,L,xx,yy,zz)))/2};
StageEndIdx = 1:length(StageEnd)-1;
end

function [PhiStack,CellIdx,StageBreak,StageEndIdx] = LoadInitState(t0,DirName,StageEnd)
[tSave, MatIdx, FileList] = GetTime(DirName);
LoadMatIdx = MatIdx(abs(tSave-t0)<1e-4);% floating-point problem
% load phi
if isempty(LoadMatIdx)
  error('.mat at %0.1f do not exist',t0)
end
FileName = sprintf('%s/%s',DirName,FileList{LoadMatIdx});
load(FileName,'phi_save','CellIdx','StageIdx')
% fprintf('load %s\n',FileName);
StageBreak = length(StageIdx) == 1;
PhiStack = phi_save;
StageEndIdx = StageIdx(1):length(StageEnd)-1;
end

function AnalResult = InitAnalVar(SimName,AnalName,LoadSwitch)
FileName = sprintf('%1$s/%1$s_%2$s.mat',SimName,AnalName);
if exist(FileName,'file')
  if LoadSwitch
    load(FileName,AnalName);
    AnalResult = eval(AnalName);
  else
    error('%s file already exist !',AnalName)
  end
else
  AnalResult = {};
end
end

function ParaSave(StrF,StrStack,phi_save,CellIdx,SimParameters,StageIdx)
unsaved = 1;
n = 0;
MatName = sprintf(StrF,StrStack{:});
while unsaved
  try
    save(MatName,'phi_save','CellIdx','SimParameters','StageIdx','-v7.3')
    unsaved = 0;
  catch ME
		disp(ME.message)
    unsaved = 1;
    n = n+1;
  end
  if n>10
    warning('%s unsaved\n',MatName)
    break
  end
end
end

function EggshellChecking(PhiStack,PhiShell,th,N)
PhaseSum = PhiStack{1};
for n = 2:length(PhiStack)
  PhaseSum = PhaseSum + PhiStack{n};
end
BoxCell = regionprops(PhaseSum>th,'BoundingBox');
BoxShell = regionprops(1-PhiShell>0.5,'BoundingBox');
if any(abs(BoxShell.BoundingBox(4:6) - BoxCell.BoundingBox(4:6))>N/2)
  warning('eggshell-cell mismatching')
end
end

function [StageBreak, SimBreak, IntervalLength, StageEndInfo, SimEndInfo, countQSS] = IsEnd(StageEndDef, DefSwitch, DefIdx, PhiStack, v, countQSS, t, t0, SimEndTime, AnalStep, h)
StageBreak = 0;
SimBreak = 0;
IntervalLength = AnalStep;
StageEndInfo = '';
SimEndInfo = '';
% SimEnd
if any(cellfun(@(x)any(isnan(x),'all'),PhiStack))% unstable
  SimEndInfo = 'Unstable';
  SimBreak = 1;
end
if any(cellfun(@(x)all(x<0.5,'all'),PhiStack))%CellLoss
  SimEndInfo = 'CellLoss';
  SimBreak = 1;
end
if t >= SimEndTime%SimTimeout
  SimEndInfo = 'Timeout';
  SimBreak = 1;
end
if SimBreak% terminate the simulation
  fprintf('%0.1f %s\n',t,SimEndInfo);
  return
end
% StageEnd
vRMS = rms(cell2mat(cellfun(@(x)vecnorm(x(max(end-2,1):end,:)',2),v,'UniformOutput',0)));
if DefSwitch(1)%Duration
  IntervalLength = ceil(min(StageEndDef{DefIdx(1),2}-t+t0, AnalStep*h)/h);
  if t-t0>=StageEndDef{DefIdx(1),2}
    StageEndInfo = 'Duration';
    StageBreak = 1;
  end
end
if DefSwitch(2)%SS
  if vRMS(end) <= StageEndDef{DefIdx(2),2}
    StageEndInfo = 'SS';
    StageBreak = 1;
  end
end
if DefSwitch(3)%QSS
	countQSS = countQSS+any(islocalmin(vRMS));
  if countQSS==StageEndDef{DefIdx(3),2}
    StageEndInfo = 'QSS';
    StageBreak = 1;
%     t = t-AnalStep*h;
  end
end
% if StageBreak
%   fprintf('%0.1f reach %s\n',t,StageEndInfo);
% end
end
