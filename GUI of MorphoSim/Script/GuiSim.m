function GuiSim(CellName,sigmaIn,Mask,dt,TimeNum,SavingInterval,SimName,EnableEggshell,EnableGPU)
tspan = [0,TimeNum];
CellParameters = saveIC(SimName,Mask,CellName);
%% enable eggshell
LongAxis = 28.284630830538536;
b = 0.665;
c = 0.455;
d = 0.361;
AxisLength = LongAxis.*[1;b;c;d];
%% enable GPU
GPUStr = {};
if EnableGPU
  GPUStr = {'GPU'};
end
%% parameter
SimParameters.tau = 2.62;
SimParameters.M = 8;
SimParameters.gamma = 0.25;
SimParameters.g = 1.6;
SimParameters.g_shell = 16;
SimParameters.c = 2;
SimParameters.sigma = sigmaIn;
% load(sprintf('%s/%s_CellParameters.mat',SimName,SimName),'CellParameters')
%% stage end
StageEnd(1).DividingCell = {};
StageEnd(1).SimParaModify = [];
StageEnd(1).StageEndDef = {'SS',1e-4};
StageEnd(2).DividingCell = {};
StageEnd(2).SimParaModify = [];
StageEnd(2).StageEndDef = {'Duration',range(tspan)};
%% simulation
dx = 0.5;
% dt = 1;
% s = 5.3;% stabilization parameter
% dt = 1.5;
s = 12;% stabilization parameter
L = [60,40,30];% domain size
N = L/dx;% grid size
a = floor(10/dt)*dt;%analysis interval
main(SimName,L,N,AxisLength,dt,SimParameters,CellParameters,StageEnd,tspan,...
  'LoadDir',SimName,...% load IC from folder named SimName
  'SavePhi',floor(SavingInterval/a),...% nearest integers divisible by analysis interval
  'StabilizationPara',s,...% stabilization parameter
  'EnforceStage',1, ...% avoid bug when load simulation.
  'EnableEggshell',EnableEggshell,...
  GPUStr{:});% enable GPU acceleration
% MorphoSim(SimName,L,N,AxisLength,dt,SimParameters,CellParameters,StageEnd,tspan,...