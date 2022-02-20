function [time, idx, FileList] = GetTime(DirName)
if nargin == 0
  DirName = '.';
end
if isfolder(DirName)
	FileNameChar = ls(DirName);
else
	time = [];
  warning('dir %s do not exist\n',DirName)
	return
end
if ispc
  FileList = string(FileNameChar);
elseif isunix
  FileList = regexp(FileNameChar,'\s','split')';
  FileList(cellfun(@isempty,FileList)) = [];
end
time = cell2mat(cellfun(@StrFun,FileList,'UniformOutput',0));
IsIC = cellfun(@(x)contains(x,'IC'),FileList);
IsMat = cellfun(@(x)contains(x,'.mat'),FileList);
IsPhi = IsMat&~(isnan(time)|IsIC);
[time,I] = sort(time(IsPhi));
idx = find(IsPhi);
idx = idx(I);
end

function out = StrFun(in)
temp = in(find(in=='_',1,'last')+1:find(in=='.')-1);
temp(temp == 'd') = '.';
out = str2double(temp);
end
