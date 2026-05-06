function dat_out = log_data(dat_in)
dat_in_sign = (dat_in);
dat_out = dat_in_sign.*log10(1+abs(dat_in));
end