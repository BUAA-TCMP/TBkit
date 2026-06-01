function Omega_ab_mu = BerryCurvature_k(Ham, tensor_index, kpoint, mu_list, options)
arguments
    Ham TBkit
    tensor_index (1,2) double
    kpoint (1,3) double
    mu_list double
    options.T = 50 % Kelvin
    options.eps = 1e-4
end
Nbands = Ham.Basis_num;
a = tensor_index(1);
b = tensor_index(2);
%%
[~, EIG_ki, ~, VEC_ki] = Ham.fft(kpoint);

dEnm = repmat(EIG_ki, 1, Nbands) - repmat(EIG_ki', Nbands, 1);
inv_dEnm = zeros(Nbands, Nbands);
is_degenerated = abs(dEnm) < options.eps;
inv_dEnm(~is_degenerated) = 1./dEnm(~is_degenerated);
%%
Omega_ab = zeros([Nbands,1]);

for n = 1:Nbands
    for m = 1:Nbands
        Omega_ab(n) = Omega_ab(n) - 2*imag(VEC_ki(n,m,a) * VEC_ki(m,n,b)) * inv_dEnm(n,m)^2;
    end
end

Nmu = length(mu_list);
E_minus_mu = repmat(EIG_ki, 1, Nmu) - repmat(mu_list, Nbands, 1);
f0 = Fermi_0(E_minus_mu, options.T);

Omega_ab_mu = tensorprod(Omega_ab, f0, 1, 1);
end