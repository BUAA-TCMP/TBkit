function H_hr = kron(A,B)
%KRON Kronecker tensor product for HR objects
%   Performs Kronecker product between HR objects or between HR and matrices.
%
%   Syntax:
%       H_hr = kron(A,B)
%
%   Inputs:
%       A - First operand (HR object or matrix)
%       B - Second operand (HR object or matrix)
%
%   Outputs:
%       H_hr - Resulting HR object from Kronecker product
%
%   Example:
%       H_kron = kron(H1, H2); % Kronecker product of two HR objects
%       H_kron = kron(H1, eye(2)); % Kronecker product with identity matrix
if isa(A,'HR') && isa(B,'HR')
H_hr1 = A;
H_hr2 = B;
H_hr =  HR(H_hr1.WAN_NUM * H_hr2.WAN_NUM,...
unique([H_hr1.vectorL;H_hr2.vectorL],'rows'));
for i = 1:H_hr.NRPTS
vector = H_hr.vectorL(i,:);
[~,seq1]=ismember(vector,H_hr1.vectorL,'rows');
[~,seq2]=ismember(vector,H_hr2.vectorL,'rows');
if seq1 ~= 0 && seq2 ~=0
amp = kron(H_hr1.HnumL(:,:,seq1) ,...
H_hr2.HnumL(:,:,seq2));
amp_sym = kron(H_hr1.HcoeL(:,:,seq1) ,...
H_hr2.HcoeL(:,:,seq2));
H_hr = H_hr.set_hop_mat(amp,vector,'set');
H_hr = H_hr.set_hop_mat(amp_sym,vector,'sym');
end
end
elseif isa(A,'HR') && ~isa(B,'HR')
H_hr1 = A;
H_hr = H_hr1;
if isa(B,'sym')
for i = 1:H_hr.NRPTS
H_hr.HcoeL(:,:,i) = kron(H_hr.HcoeL(:,:,i),B);
end
elseif isa(B,'numeric')
for i = 1:H_hr.NRPTS
H_hr.HnumL(:,:,i) = kron(H_hr.HnumL(:,:,i),B);
end
else
end
elseif ~isa(A,'HR') && isa(B,'HR')
H_hr1 = B;
H_hr = H_hr1;
if isa(A,'sym')
for i = 1:H_hr.NRPTS
H_hr.HcoeL(:,:,i) = kron( A , H_hr.HcoeL(:,:,i));
end
elseif isa(A,'numeric')
for i = 1:H_hr.NRPTS
H_hr.HnumL(:,:,i) = kron( A , H_hr.HnumL(:,:,i));
end
else
end
end
end
