function H_hr = mrdivide(A,B)
% MRDIVIDE Overloaded right matrix division for HR objects
%
%   H_HR = MRDIVIDE(A,B) implements right matrix division (A/B) for HR objects
%
%   Inputs:
%       A - HR object or numeric matrix
%       B - HR object or numeric matrix
%   Output:
%       H_hr - Resulting HR object after division
%
%   Notes:
%       - Handles three cases: HR/HR, HR/numeric, numeric/HR
%       - Checks WAN_NUM compatibility for HR/HR case
%       - Supports both numeric and symbolic operations

if isa(A,'HR') && isa(B,'HR')
H_hr1 = A;
H_hr2 = B;
if H_hr1.WAN_NUM ~= H_hr2.WAN_NUM
error('WAN_NUM different');
end
H_hr =  HR(H_hr2.WAN_NUM,...
unique([H_hr1.vectorL;H_hr2.vectorL],'rows'));
for i = 1:H_hr.NRPTS
vector = H_hr.vectorL(i,:);
[~,seq1]=ismember(vector,H_hr1.vectorL,'rows');
[~,seq2]=ismember(vector,H_hr2.vectorL,'rows');
if seq1 ~= 0 && seq2 ~=0
amp = H_hr1.HnumL(:,:,seq1) / H_hr2.HnumL(:,:,seq2);
amp_sym = H_hr1.HcoeL(:,:,seq1) / H_hr2.HcoeL(:,:,seq2);
H_hr = H_hr.set_hop_mat(amp,vector,'set');
H_hr = H_hr.set_hop_mat(amp_sym,vector,'sym');
end
end
elseif isa(A,'HR') && ~isa(B,'HR')
H_hr = A;
if A.Type == 'list'
H_hr.HcoeL = H_hr.HcoeL/B;
H_hr.HnumL = H_hr.HnumL/B;
else
if isa(B,'sym')
for i = 1:H_hr.NRPTS
H_hr.HcoeL(:,:,i) = H_hr.HcoeL(:,:,i) / B;
end
elseif isa(B,'numeric')
for i = 1:H_hr.NRPTS
H_hr.HnumL(:,:,i) = H_hr.HnumL(:,:,i) / B;
end
else
end
end
elseif ~isa(A,'HR') && isa(B,'HR')
H_hr1 = B;
H_hr = H_hr1;
H_hr = H_hr.line_000_gen();
i = H_hr.Line_000;
if H_hr.WAN_NUM ~= length(B)
error('WAN_NUM different');
end
if isa(A,'sym')
H_hr.HcoeL(:,:,i) = A / H_hr.HcoeL(:,:,i) ;
elseif isa(A,'numeric')
H_hr.HnumL(:,:,i) = A / H_hr.HnumL(:,:,i) ;
else
end
end
end
