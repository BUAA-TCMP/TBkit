function H_hk = hermitize(H_hk)
%HERMITIZE Enforce Hermiticity on Hamiltonian
%
% Syntax:
%   H_hk = hermitize(H_hk)
%
% Input:
%   H_hk - HK Hamiltonian object
%
% Output:
%   H_hk - Hermitized Hamiltonian
%
% Description:
%   Ensures Hamiltonian is Hermitian by averaging with its conjugate
%   transpose. Handles both:
%   - Numeric coefficients (HnumL)
%   - Symbolic coefficients (HcoeL)
%
% Algorithm:
%   H_herm = (H + H')/2
%
% Note:
%   Preserves original values if already Hermitian
%   Automatically skips empty coefficient sets
%
% Example:
%   H_herm = Hk_obj.hermitize();
if isequal(zeros(size(H_hk.HnumL)),H_hk.HnumL)
    num_label = false;
else
    num_label = true;
end
if isequal(sym(zeros(size(H_hk.HcoeL))),H_hk.HcoeL)
    coe_label = false;
else
    coe_label = true;
end
H_hk_bk = H_hk';
if coe_label
    H_hk = (H_hk + H_hk_bk)/2;
end
if num_label
    H_hk.HnumL = (H_hk_bk.HnumL + H_hk.HnumL )/2;
end
end
