function result = TRAHC_int_ik_2D(Ham, kpoint, mu_list, T_list, selectbands, varargin)
% TRAHC_int_ik_2D
%
% Single-k intrinsic third-order anomalous Hall integrand in 2D.
%
% This is the MATLAB analogue of the Fortran routine TRAHC_int_k.
% It follows the same architecture as ForthHall_ik_2D_sym_sumrule:
%
%   1. call Tabc_k_fast_2D(Ham, kpoint, varargin{:})
%   2. extract band energies, diagonal velocities and T^{abc}
%   3. loop over chemical potential and Fermi-smearing scale
%   4. accumulate
%
% Formula implemented:
%
%   chi_{a b g e} = f0'(epsilon_n - mu)
%                 * [ v_n^a T_n^{b g e} - v_n^b T_n^{a g e} ]
%
% The Fortran-exported components are
%
%   tensor_out(:,:,1,1) = chi_yxxx = chi_{2111}
%   tensor_out(:,:,1,2) = chi_yxxy = chi_{2112}
%   tensor_out(:,:,1,3) = chi_yxyx = chi_{2121}
%   tensor_out(:,:,1,4) = chi_yxyy = chi_{2122}
%
% Inputs
% ------
% Ham         : Hamiltonian object
% kpoint      : [kx, ky, kz]
% mu_list     : chemical potential list, same role as Ef_list in Fortran
% T_list      : Fermi-smearing / temperature list, same role as Eta_array in Fortran
% selectbands : selected band indices; if empty, use all bands
% varargin    : name-value pairs passed to Tabc_k_fast_2D, e.g. 'eta',1e-5
%
% Outputs
% -------
% result.tensor_out      : [nmu, nT, 1, 4], exactly matching Fortran ordering
% result.tensor_raw      : [nmu, nT, 2,2,2,2], full unsymmetrized chi_{a b g e}
% result.tensor_yx_raw   : [nmu, nT, 4], [yxxx, yxxy, yxyx, yxyy]
% result.tensor_yx_sym_ge: [nmu, nT, 3], [yxxx, yx(xy), yxyy] with g/e symmetrized
%
% Notes
% -----
% - The first two indices a,b are kept explicitly antisymmetric by the formula.
% - The last two indices g,e are NOT symmetrized in tensor_out, to reproduce Fortran.
%   If the two external electric fields are physically indistinguishable, use
%   result.tensor_yx_sym_ge(:, :, 2) = 0.5*(chi_yxxy + chi_yxyx).
% - This function returns the local k-resolved integrand only. Brillouin-zone
%   weights, prefactors, spin degeneracy factors, and MPI/k-loop accumulation
%   should be applied outside, just as in the Fortran driver.

if nargin < 5 || isempty(selectbands)
    selectbands = 1:Ham.Basis_num;
end

mu_list = mu_list(:);
T_list  = T_list(:).';

if any(T_list <= 0)
    error('TRAHC_int_ik_2D:BadTemperature', ...
          'All entries of T_list must be positive.');
end

% ------------------------------------------------------------
% Evaluate single-k band quantities from Tabc_k_fast_2D
% ------------------------------------------------------------
out = Tabc_k_fast_2D(Ham, kpoint, varargin{:});

sel  = selectbands(:).';
nsel = numel(sel);
nmu  = numel(mu_list);
nT   = numel(T_list);

% band energies
E = reshape(out.eigvals(sel), [], 1);      % [nsel,1]

% diagonal velocities v_a^{nn}, a = x,y
Vdiag = zeros(2, nsel);
for a = 1:2
    Va_full = reshape(out.vel(a,:,:), size(out.vel,2), size(out.vel,3));
    Vdiag(a,:) = real(diag(Va_full(sel, sel))).';
end

% T_tensor(a,b,c,n)
Tten = out.T_tensor(:,:,:,sel);            % [2,2,2,nsel]

% ------------------------------------------------------------
% Allocate outputs
% ------------------------------------------------------------
tensor_out = zeros(nmu, nT, 1, 4);
tensor_raw = zeros(nmu, nT, 2, 2, 2, 2);

% ------------------------------------------------------------
% Main loops: mu, T/Fermi-smearing, tensor indices, bands
% ------------------------------------------------------------
for iT = 1:nT
    Tnow = T_list(iT);

    for imu = 1:nmu
        mu = mu_list(imu);
        f0p = Fermi_1_fortran_local(E - mu, Tnow);  % [nsel,1]
        f0p_row = f0p.';                            % [1,nsel]

        for a = 1:2
            va = Vdiag(a,:);

            for b = 1:2
                vb = Vdiag(b,:);

                for g = 1:2
                    for e = 1:2
                        T_bge = reshape(Tten(b,g,e,:), 1, nsel);
                        T_age = reshape(Tten(a,g,e,:), 1, nsel);

                        chi_abge = sum(f0p_row .* (va .* T_bge - vb .* T_age));
                        tensor_raw(imu, iT, a, b, g, e) = chi_abge;
                    end
                end
            end
        end

        % Fortran component order: [yxxx, yxxy, yxyx, yxyy]
        tensor_out(imu, iT, 1, 1) = tensor_raw(imu, iT, 2, 1, 1, 1);
        tensor_out(imu, iT, 1, 2) = tensor_raw(imu, iT, 2, 1, 1, 2);
        tensor_out(imu, iT, 1, 3) = tensor_raw(imu, iT, 2, 1, 2, 1);
        tensor_out(imu, iT, 1, 4) = tensor_raw(imu, iT, 2, 1, 2, 2);
    end
end

% convenient views
tensor_yx_raw = zeros(nmu, nT, 4);
tensor_yx_raw(:,:,1) = tensor_out(:,:,1,1);    % yxxx
tensor_yx_raw(:,:,2) = tensor_out(:,:,1,2);    % yxxy
tensor_yx_raw(:,:,3) = tensor_out(:,:,1,3);    % yxyx
tensor_yx_raw(:,:,4) = tensor_out(:,:,1,4);    % yxyy

tensor_yx_sym_ge = zeros(nmu, nT, 3);
tensor_yx_sym_ge(:,:,1) = tensor_yx_raw(:,:,1);
tensor_yx_sym_ge(:,:,2) = 0.5 * (tensor_yx_raw(:,:,2) + tensor_yx_raw(:,:,3));
tensor_yx_sym_ge(:,:,3) = tensor_yx_raw(:,:,4);

% ------------------------------------------------------------
% Pack output
% ------------------------------------------------------------
result = struct();
result.kpoint      = kpoint(:).';
result.mu_list     = mu_list;
result.T_list      = T_list;
result.selectbands = sel;

result.tensor_out       = tensor_out;
result.tensor_raw       = tensor_raw;
result.tensor_yx_raw    = tensor_yx_raw;
result.tensor_yx_sym_ge = tensor_yx_sym_ge;

result.component_names_raw    = {'chi_yxxx','chi_yxxy','chi_yxyx','chi_yxyy'};
result.component_names_sym_ge = {'chi_yxxx','chi_yxxy_sym','chi_yxyy'};

% Also expose ingredients for debugging/benchmarking against Fortran
result.eigvals = E;
result.vdiag   = Vdiag;

end

% ============================================================
% local helper: exact Fortran-style f0'(epsilon - mu)
% ============================================================
function f1 = Fermi_1_fortran_local(x, T)
% For f = 1/(exp(x/T)+1), return df/dx.
% Same cutoff convention as the Fortran code: |x/T| > 50 -> 0.

f1 = zeros(size(x));
u = x ./ T;
mask = abs(u) <= 50.0;

if any(mask(:))
    F = 1.0 ./ (exp(u(mask)) + 1.0);
    f1(mask) = -F .* (1.0 - F) ./ T;
end

end
