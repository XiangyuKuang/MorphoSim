function CellParameters = saveIC(SimName,Mask,CellName)
if exist(SimName,'dir') ~= 7% make dir named by SimName
   mkdir(SimName)
end
% save IC
phi_save = cellfun(@(x)Mask2PF(double(x)),Mask,'UniformOutput',0);
CellIdx = (1:length(phi_save))';
StageIdx = 1;% do not delete this line
save(sprintf('%s/%s_0.mat',SimName,SimName),'phi_save','CellIdx','StageIdx')
% save CellPara
CellParameters.Name = CellName;
CellParameters.Volume = cellfun(@(x)sum(x,'all')*0.5^3,phi_save);
% save(sprintf('%s/%s_CellParameters.mat',SimName,SimName),'CellParameters')