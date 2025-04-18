function H_hr = rmhop(H_hr,vectorL,options)
% RMHOP Remove hopping terms from HR object
%
%   H_HR = RMHOP(H_HR,VECTORL,OPTIONS) removes specified hopping terms
%
%   Inputs:
%       H_hr - HR object to modify
%       vectorL - Vectors to remove (can be indices or symbolic patterns)
%       options.enforce_list - Keep in list format [default: false]
%   Output:
%       H_hr - Modified HR object
%
%   Notes:
%       - Can remove by vector indices or symbolic pattern
%       - Optionally maintains list format
%       - Preserves remaining hopping terms
arguments
H_hr HR
vectorL = [1 1];
options.enforce_list  = false;
end
if ~strcmp(H_hr.type, 'list')
H_hr = H_hr.rewrite();
end
if options.enforce_list
giveback  = false;
else
giveback  = true;
end
vectorList = double(H_hr.vectorL);
if isa(vectorL,'double')
switch size(vectorL,2)
case 2
[seq] = ~(ismember(H_hr.vectorL(:,H_hr.Dim+1:H_hr.Dim+2),(vectorL),"rows"));
case 3
[seq] = ~(ismember(H_hr.vectorL(:,1:3),(vectorL),"rows"));
case 5
[seq] = ~(ismember(H_hr.vectorL(:,1:5),(vectorL),"rows"));
end
elseif isa(vectorL,'sym')
[seq] = find(park.strcontain(string(H_hr.HcoeL),string(vectorL)));
end
if H_hr.num
H_hr.HnumL=H_hr.HnumL(seq,:);
end
if H_hr.coe
H_hr.HcoeL=H_hr.HcoeL(seq,:);
end
H_hr.vectorL=H_hr.vectorL(seq,:);
if giveback
H_hr = H_hr.rewind();
else
H_hr = H_hr;
end
end
