function [Omega_ab,VEC_ki] = BerryCurvature_nk(Ham, tensor_index, kpoint, options)
arguments
    Ham TBkit
    tensor_index (1,2) double
    kpoint (1,3) double
    options.eps = 1e-4
    options.rotate_cart = [];
end
Nbands = Ham.Basis_num;
a = tensor_index(1);
b = tensor_index(2);
%%
if isempty(options.rotate_cart)
    [~,EIG_ki,~,VEC_ki] = Ham.fft(kpoint);
else
    [~,EIG_ki,~,VEC_ki] = Ham.fft(kpoint,options.rotate_cart);
end


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
end