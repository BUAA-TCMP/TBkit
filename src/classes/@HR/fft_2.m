function [W,D,dH_dk_R,dH_dk_dk_R,Hk] = fft_2(H_hr, klist_cart)
arguments
    H_hr HR
    klist_cart (:,1) double   % Cartesian k-vector
end

% =========================================================
% Fast Bloch Hamiltonian with orbital-position correction
% fractional orbital coordinates are assumed:
%
% H_ij(k) = sum_R H_ij(R) exp[i k·(R + tau_j - tau_i)]
%
% where:
%   R_cart   = vectorL * Rm
%   tau_cart = orbL    * Rm
%
% Speed strategy:
%   exp[i k·(R + dtau_ij)] = exp(i k·R) * exp(i k·dtau_ij)
% so keep the fast tensorprod over R only, then multiply by ij-phase.
% =========================================================

persistent cache

Dim = H_hr.Dim;
k = klist_cart(:);
assert(length(k) == Dim, 'klist_cart dimension mismatch.');

HnumList = H_hr.HnumL;   % norb x norb x nR
[norb, norb2, nR] = size(HnumList);
assert(norb == norb2, 'HnumL must be norb x norb x nR.');

% ---------------------------------------------------------
% cache rebuild
% ---------------------------------------------------------
need_rebuild = true;
if ~isempty(cache)
    if isfield(cache,'norb') && isfield(cache,'nR') && isfield(cache,'Dim')
        if cache.norb == norb && cache.nR == nR && cache.Dim == Dim
            need_rebuild = false;
        end
    end
end

if need_rebuild
    % ----- R vectors in Cartesian -----
    vectorList = double(H_hr.vectorL(:,1:Dim));   % nR x Dim (fractional lattice coords)
    R_cart = vectorList * H_hr.Rm;                % nR x Dim

    % ----- orbital positions in Cartesian -----
    % USER SAID: orbL is fractional coordinates
    tau_frac = double(H_hr.orbL(:,1:Dim));        % norb x Dim
    tau_cart = tau_frac * H_hr.Rm;                % norb x Dim

    % ----- dtau(i,j,a) = tau_j(a) - tau_i(a) -----
    dtau = zeros(norb, norb, Dim);
    for a = 1:Dim
        dtau(:,:,a) = tau_cart(:,a).' - tau_cart(:,a);
    end

    % ----- pack cache -----
    cache.R_cart = R_cart;
    cache.dtau   = dtau;
    cache.norb   = norb;
    cache.nR     = nR;
    cache.Dim    = Dim;
end

R_cart = cache.R_cart;   % nR x Dim
dtau   = cache.dtau;     % norb x norb x Dim

% ---------------------------------------------------------
% phase over R only
% ---------------------------------------------------------
phase_R = exp(1i * (R_cart * k));     % nR x 1

% S0 = sum_R H(R) e^{ik·R}
S0 = tensorprod(HnumList, phase_R, 3, 1);   % norb x norb

% S1_a = sum_R R_a H(R) e^{ik·R}
S1 = zeros(norb, norb, Dim);
for a = 1:Dim
    S1(:,:,a) = tensorprod(HnumList, phase_R .* R_cart(:,a), 3, 1);
end

% S2_ab = sum_R R_a R_b H(R) e^{ik·R}
S2 = zeros(norb, norb, Dim, Dim);
for a = 1:Dim
    Ra = R_cart(:,a);
    for b = 1:Dim
        S2(:,:,a,b) = tensorprod(HnumList, phase_R .* Ra .* R_cart(:,b), 3, 1);
    end
end

% ---------------------------------------------------------
% orbital phase P_ij = exp[i k·(tau_j - tau_i)]
% ---------------------------------------------------------
phase_tau_arg = zeros(norb, norb);
for a = 1:Dim
    phase_tau_arg = phase_tau_arg + k(a) * dtau(:,:,a);
end
P_tau = exp(1i * phase_tau_arg);   % norb x norb

% ---------------------------------------------------------
% H(k)
% ---------------------------------------------------------
Hk = P_tau .* S0;

% ---------------------------------------------------------
% first derivative
% dH_a = i P_tau [ S1_a + dtau_a S0 ]
% ---------------------------------------------------------
dH_dk_R = zeros(norb, norb, Dim);
for a = 1:Dim
    dH_dk_R(:,:,a) = 1i * P_tau .* ( S1(:,:,a) + dtau(:,:,a) .* S0 );
end

% ---------------------------------------------------------
% second derivative
% d2H_ab = - P_tau [ S2_ab + dtau_a S1_b + dtau_b S1_a + dtau_a dtau_b S0 ]
% ---------------------------------------------------------
dH_dk_dk_R = zeros(norb, norb, Dim, Dim);
for a = 1:Dim
    dta = dtau(:,:,a);
    for b = 1:Dim
        dtb = dtau(:,:,b);
        dH_dk_dk_R(:,:,a,b) = - P_tau .* ( ...
            S2(:,:,a,b) ...
            + dta .* S1(:,:,b) ...
            + dtb .* S1(:,:,a) ...
            + dta .* dtb .* S0 );
    end
end

% ---------------------------------------------------------
% Hermitian symmetrization
% ---------------------------------------------------------
Hk = (Hk + Hk') / 2;
for a = 1:Dim
    dH_dk_R(:,:,a) = (dH_dk_R(:,:,a) + dH_dk_R(:,:,a)') / 2;
    for b = 1:Dim
        dH_dk_dk_R(:,:,a,b) = (dH_dk_dk_R(:,:,a,b) + dH_dk_dk_R(:,:,a,b)') / 2;
    end
end

% ---------------------------------------------------------
% eig
% ---------------------------------------------------------
[W,D] = eig(Hk, 'vector');

end

% function [W,D,dH_dk_R,dH_dk_dk_R] = fft_2(H_hr, klist_cart)
% arguments
%     H_hr HR
%     klist_cart % cart
% end
% 
% % try
% %     [W, D, dH_dk_R] = fft_mex(H_hr.HnumL, double(H_hr.vectorL(:,1:H_hr.Dim))* H_hr.Rm, klist);
% % catch
% %     HnumList = H_hr.HnumL ;
% %     vectorList = double(H_hr.vectorL(:,1:H_hr.Dim)) ;
% %     vectorList_R = vectorList * H_hr.Rm;
% %     FactorList = exp(1i*vectorList_R*klist.');
% %     Hout = tensorprod(HnumList, FactorList, 3, 1);
% %     [W,D]= eig((Hout+Hout')/2,'vector');
% %     for i = 1:3
% %         dH_dk_R(:,:,i) = 1i * tensorprod(HnumList , FactorList.*vectorList_R(:,i), 3, 1); % only for one kpoint?
% %     end
% % end
% % % return;
% persistent HnumList vectorList_R 
% if isempty(HnumList)
%     HnumList = H_hr.HnumL;
%     vectorList = double(H_hr.vectorL(:,1:H_hr.Dim));  % M×3
%     vectorList_R = vectorList * H_hr.Rm;  % M×3 real-space R
% end
% 
% 
% % 计算因子：exp(i * R · k)
% % orbL
% phase = exp(1i * (vectorList_R * klist_cart(:)));  % M×1
% Hk = tensorprod(HnumList, phase, 3, 1);
% Hk = (Hk + Hk') / 2;  % Hermitian 修正
% [W,D] = eig(Hk, 'vector');
% for i = 1:3
%     dH_dk_R(:,:,i) = 1i * tensorprod(HnumList , phase.*vectorList_R(:,i), 3, 1);
%     for j = 1:3
%         dH_dk_dk_R(:,:,i,j) = - tensorprod(HnumList , phase.*vectorList_R(:,i).*vectorList_R(:,j), 3, 1);
%     end
% end
% 
% 
% end
