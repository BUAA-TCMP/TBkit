function Omega_ab_mu = BCD_k(Ham, tensor_index, kpoint, mu_list, options)
arguments
    Ham TBkit
    tensor_index (1,3) double
    kpoint (1,3) double
    mu_list double
    options.T = 50 % Kelvin
    options.eps = 1e-4
end
Nbands = Ham.Basis_num;
a = tensor_index(1);
b = tensor_index(2);
c = tensor_index(3);
%%
[WAV_ki, EIG_ki, dH_dk_xyz] = Ham.fft(kpoint);

dEnm = repmat(EIG_ki, 1, Nbands) - repmat(EIG_ki', Nbands, 1);
inv_dEnm = zeros(Nbands, Nbands);
is_degenerated = abs(dEnm) < options.eps;
inv_dEnm(~is_degenerated) = 1./dEnm(~is_degenerated);
%%
VEC_ki = zeros(Nbands, Nbands, 3);
for i = 1:3
    VEC_ki(:,:,i) = WAV_ki' * dH_dk_xyz(:,:,i) * WAV_ki;
end
%%
Omega_ab = zeros([Nbands,1]);

for n = 1:Nbands
    for m = 1:Nbands
        Omega_ab(n) = Omega_ab(n) - 2*imag(VEC_ki(n,m,a) * VEC_ki(m,n,b)) * inv_dEnm(n,m)^2;
    end
end

Nmu = length(mu_list);
E_minus_mu = repmat(EIG_ki, 1, Nmu) - repmat(mu_list, Nbands, 1);
f1 = Fermi_1(E_minus_mu, options.T);

Omega_ab_mu = tensorprod(Omega_ab .* real(diag(VEC_ki(:,:,c))), f1, 1, 1);
end