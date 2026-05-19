function result = ForthHall_ik_2D_sym_sumrule(Ham, kpoint, mu_list, T, tau, selectbands, varargin)
% ForthHall_ik_2D_sym_sumrule
%
% Single-k 2D local response using sum-rule / integration-by-parts form.
% Directly calls:
%   out = Tabc_k_fast_2D(Ham, kpoint, varargin{:})
%
% Inputs
% ------
% Ham          : Hamiltonian object
% kpoint       : [kx, ky, kz]
% mu_list      : list of chemical potentials
% T            : temperature (must be > 0 for smooth Fermi derivatives)
% tau          : relaxation time
% selectbands  : selected band indices; if empty, use all bands
% varargin     : name-value pairs passed to Tabc_k_fast_2D
%
% Output
% ------
% result struct with fields:
%   .chi_tau1   : [nmu, 6]
%   .chi_tau3   : [nmu, 6]
%   .chi_total  : [nmu, 6]
%   .sym_tau1   : [nmu, 2,2,2,2,2]
%   .sym_tau3   : [nmu, 2,2,2,2,2]
%   .sym_total  : [nmu, 2,2,2,2,2]
%   .vec_tau1   : [nmu, 10]
%   .vec_tau3   : [nmu, 10]
%   .vec_total  : [nmu, 10]

if nargin < 6 || isempty(selectbands)
    selectbands = 1:Ham.Basis_num;
end

if T <= 0
    error('This sum-rule implementation assumes finite temperature T > 0.');
end

% ------------------------------------------------------------
% Evaluate single-point band quantities from Tabc_k_fast_2D
% ------------------------------------------------------------
out = Tabc_k_fast_2D(Ham, kpoint, varargin{:});

sel  = selectbands(:).';
nsel = numel(sel);
nmu  = numel(mu_list);

% eigvals
E = out.eigvals(sel);   % [nsel,1]

% diagonal velocities v_a^{nn}
Vdiag = zeros(2, nsel);
for a = 1:2
    Va_full = reshape(out.vel(a,:,:), size(out.vel,2), size(out.vel,3));
    dv_full = diag(Va_full);
    Vdiag(a,:) = real(dv_full(sel)).';
end

% dvel_diag(a,b,n)
dvel = out.dvel_diag(:,:,sel);          % [2,2,nsel]

% T_tensor(c,d,e,n)
Tten = out.T_tensor(:,:,:,sel);         % [2,2,2,nsel]

% BCP / G_{de}(n)
Gten = out.BCP(:,:,sel);                % [2,2,nsel]

% dOmega(a,b,c,n)
dOm  = out.dOmega(:,:,:,sel);           % [2,2,2,nsel]

% ------------------------------------------------------------
% Preallocate outputs
% ------------------------------------------------------------
sym_tau1  = zeros(nmu, 2,2,2,2,2);
sym_tau3  = zeros(nmu, 2,2,2,2,2);
sym_total = zeros(nmu, 2,2,2,2,2);

vec_tau1  = zeros(nmu, 10);
vec_tau3  = zeros(nmu, 10);
vec_total = zeros(nmu, 10);

chi_tau1  = zeros(nmu, 6);
chi_tau3  = zeros(nmu, 6);
chi_total = zeros(nmu, 6);

% ordered independent bcde sectors for 2D
ordered_bcdes = [
    1 1 1 1
    1 1 1 2
    1 1 2 2
    1 2 2 2
    2 2 2 2
];

% ------------------------------------------------------------
% Loop over mu
% ------------------------------------------------------------
for imu = 1:nmu
    mu = mu_list(imu);

    Emu  = E - mu;                 % [nsel,1]
    f0p  = Fermi_1_local(Emu, T);  % [nsel,1]
    f0pp = Fermi_2_local(Emu, T);  % [nsel,1]

    % unsymmetrized local integrands
    raw_tau1 = zeros(2,2,2,2,2);
    raw_tau3 = zeros(2,2,2,2,2);

    for a = 1:2
        va = Vdiag(a,:);   % [1,nsel]

        for b = 1:2
            vb = Vdiag(b,:);

            for c = 1:2
                dOmega_abc = reshape(dOm(a,b,c,:), 1, nsel);

                for d = 1:2
                    vd = Vdiag(d,:);

                    for e = 1:2
                        ve = Vdiag(e,:);

                        % ---- pieces used repeatedly ----
                        T_cde = reshape(Tten(c,d,e,:), 1, nsel);
                        T_bcd = reshape(Tten(b,c,d,:), 1, nsel);
                        T_acd = reshape(Tten(a,c,d,:), 1, nsel);

                        dvel_ba = reshape(dvel(b,a,:), 1, nsel);
                        dvel_ae = reshape(dvel(a,e,:), 1, nsel);
                        dvel_be = reshape(dvel(b,e,:), 1, nsel);
                        dvel_ed = reshape(dvel(e,d,:), 1, nsel);

                        G_de = reshape(Gten(d,e,:), 1, nsel);

                        % ==================================================
                        % tau1 part
                        % ==================================================
                        term_tau1 = ...
                              (-2/3) * tau * sum(va .* vb .* T_cde .* f0pp.') ...
                            + (-4/3) * tau * sum(dvel_ba .* T_cde .* f0p.') ...
                            + (   1  )* tau * sum(dvel_ae .* T_bcd .* f0p.') ...
                            - (    1 ) *tau * sum(dvel_be .* T_acd .* f0p.') ...
                            + (     1) *tau * sum(va .* ve .* T_bcd .* f0pp.') ...
                            - (     1) *tau * sum(vb .* ve .* T_acd .* f0pp.') ...
                            - (1/2) * tau * sum(dOmega_abc .* G_de .* f0p.');

                        raw_tau1(a,b,c,d,e) = term_tau1;

                        % ==================================================
                        % tau3 part
                        % ==================================================
                        term_tau3 = ...
                              tau^3 * sum(dOmega_abc .* dvel_ed .* f0p.') ...
                            + tau^3 * sum(dOmega_abc .* vd .* ve .* f0pp.');

                        raw_tau3(a,b,c,d,e) = term_tau3;
                    end
                end
            end
        end
    end

    % total
    raw_total = raw_tau1 + raw_tau3;

    % --------------------------------------------------------
    % Symmetrize over (b,c,d,e)
    % --------------------------------------------------------
    sym1 = zeros(size(raw_tau1));
    sym3 = zeros(size(raw_tau3));
    symt = zeros(size(raw_total));

    for a = 1:2
        for iord = 1:size(ordered_bcdes,1)
            bcde = ordered_bcdes(iord,:);
            plist = unique(perms(bcde), 'rows');

            acc1 = 0;
            acc3 = 0;
            acct = 0;

            for ip = 1:size(plist,1)
                b = plist(ip,1);
                c = plist(ip,2);
                d = plist(ip,3);
                e = plist(ip,4);

                acc1 = acc1 + raw_tau1(a,b,c,d,e);
                acc3 = acc3 + raw_tau3(a,b,c,d,e);
                acct = acct + raw_total(a,b,c,d,e);
            end

            acc1 = acc1 / size(plist,1);
            acc3 = acc3 / size(plist,1);
            acct = acct / size(plist,1);

            for ip = 1:size(plist,1)
                b = plist(ip,1);
                c = plist(ip,2);
                d = plist(ip,3);
                e = plist(ip,4);

                sym1(a,b,c,d,e) = acc1;
                sym3(a,b,c,d,e) = acc3;
                symt(a,b,c,d,e) = acct;
            end
        end
    end

    sym_tau1(imu,:,:,:,:,:)  = sym1;
    sym_tau3(imu,:,:,:,:,:)  = sym3;
    sym_total(imu,:,:,:,:,:) = symt;

    % --------------------------------------------------------
    % Flatten 10 independent components in the same order
    % as your old MATLAB code
    % --------------------------------------------------------
    vec1 = flatten_10_local(sym1);
    vec3 = flatten_10_local(sym3);
    vect = flatten_10_local(symt);

    vec_tau1(imu,:)  = vec1;
    vec_tau3(imu,:)  = vec3;
    vec_total(imu,:) = vect;

    % --------------------------------------------------------
    % Compose the same 6 independent chi_i as before
    % --------------------------------------------------------
    chi_tau1(imu,:)  = vec10_to_chi6_local(vec1);
    chi_tau3(imu,:)  = vec10_to_chi6_local(vec3);
    chi_total(imu,:) = vec10_to_chi6_local(vect);
end

% ------------------------------------------------------------
% pack output
% ------------------------------------------------------------
result = struct();
result.kpoint    = kpoint(:).';
result.mu_list   = mu_list(:);
result.temperature = T;
result.tau       = tau;
result.selectbands = sel;

result.chi_tau1  = chi_tau1;
result.chi_tau3  = chi_tau3;
result.chi_total = chi_total;

result.vec_tau1  = vec_tau1;
result.vec_tau3  = vec_tau3;
result.vec_total = vec_total;

result.sym_tau1  = sym_tau1;
result.sym_tau3  = sym_tau3;
result.sym_total = sym_total;

end

% ============================================================
% local helpers
% ============================================================

function f1 = Fermi_1_local(x, T)
% f'(x) for f = 1/(exp(x/T)+1)
u  = x ./ (2*T);
cu = cosh(u);
f1 = -(1 ./ (4*T)) ./ (cu.^2);
end

function f2 = Fermi_2_local(x, T)
% f''(x) for f = 1/(exp(x/T)+1)
u  = x ./ (2*T);
cu = cosh(u);
tu = tanh(u);
f2 = (1 ./ (4*T^2)) .* (1 ./ (cu.^2)) .* tu;
end

function vec = flatten_10_local(sym_raw)
% same ordering as your old code:
% a = 1,2 ; b<=c<=d<=e
vec = zeros(1,10);
idx = 1;
for a = 1:2
    for b = 1:2
        for c = b:2
            for d = c:2
                for e = d:2
                    vec(idx) = sym_raw(a,b,c,d,e);
                    idx = idx + 1;
                end
            end
        end
    end
end
end

function chi6 = vec10_to_chi6_local(vec)
% same post-processing as your old MATLAB code
chi1 = 4*vec(8)  - vec(1);
chi2 = 6*vec(9)  - 4*vec(2);
chi3 = 4*vec(10) - 6*vec(3);
chi4 = vec(6)    - 4*vec(4);
chi5 = -vec(5);
chi6 = -vec(7);
chi6 = [chi1, chi2, chi3, chi4, chi5, chi6];
end