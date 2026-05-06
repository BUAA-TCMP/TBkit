function gen_list_xyz = Magnetic_Sym_El_read(filename, BNS_number1, options)
% read generators of a magnetic space group from Mvasp2trace format file
% the output generators are all in Cartesian coordinate system 
arguments
    filename string
    BNS_number1 = 0
    options.read_translation logical = false
end

fid = fopen(filename);
line_1 = fgetl(fid);
oper_num = str2double(line_1);
if BNS_number1 == 0
    error("You must specify the space group number in order to judge the type of Bravais lattices")
end

R_abc2xyz = BCbasis2xyz(BNS_number1);
%%
gen_list_xyz = repmat(Oper(), 1, oper_num);
for i = 1:oper_num
    line_i = fgetl(fid);
    oper_i_list_form = str2num(line_i); % must use str2num here
    oper_i_isUnitary = oper_i_list_form(end)==1;

    oper_i_R = reshape(oper_i_list_form(1:9),[3,3])';
    oper_i_R = R_abc2xyz * oper_i_R * R_abc2xyz^(-1);

    if options.read_translation
        oper_i_t = oper_i_list_form(10:12);

    else
        oper_i_t = [0 0 0];
    end
    oper_i = Oper(double(oper_i_R), NaN, oper_i_t, "conjugate",~oper_i_isUnitary);
    gen_list_xyz(i) = oper_i;
end

fclose(fid);
end