function cif_struct = cif_read(filename)
% for the easy axis tag, we have assumed it is a collinear magnet
if ~contains(filename, "cif") && ~contains(filename, "mcif")
    error("The input file must be *.cif or *.mcif")
end
%% normal cif
fid = fopen(filename, 'r');
while ~feof(fid)
    line = fgetl(fid);
    
    if contains(line, '_cell_length_a')
        line = split(line);
        cif_struct.a = str2double(str_simplify(line{2}));
        continue
    end
    
    if contains(line, '_cell_length_b')
        line = split(line);
        cif_struct.b = str2double(str_simplify(line{2}));
        continue
    end

    if contains(line, '_cell_length_c')
        line = split(line);
        cif_struct.c = str2double(str_simplify(line{2}));
        continue
    end

    if contains(line, '_cell_angle_alpha')
        line = split(line);
        cif_struct.alpha = str2double(str_simplify(line{2}));
        continue
    end

    if contains(line, '_cell_angle_beta')
        line = split(line);
        cif_struct.beta = str2double(str_simplify(line{2}));
        continue
    end

    if contains(line, '_cell_angle_gamma')
        line = split(line);
        cif_struct.gamma = str2double(str_simplify(line{2}));
        continue
    end
end
frewind(fid);

Rm = abc2Rm(cif_struct.a, cif_struct.b, cif_struct.c, ...
    cif_struct.alpha, cif_struct.beta, cif_struct.gamma);

cif_struct.Rm = Rm;
R_abc2xyz = (Rm ./ vecnorm(Rm, 2, 2))';
%% symop
line_symop_start = false;
symop_xyz = Oper(eye(3));
syms x y z real;
vars = [x,y,z];

while ~feof(fid)
    line = fgetl(fid);
    
    if ~line_symop_start
        if contains(line, 'symop') && ~contains(line, 'centering')
            line_symop_start = true;
        end
    else
        if isempty(line)
            break
        elseif line(1)=="_"
            continue
        else
            line_operation = split(line);
            operation_info = split(line_operation{2},',');

            operation_info(1:3) = regexprep(operation_info(1:3), '(\d)([a-z])', '$1*$2'); % 2a -> 2*a
  
            rotation_mat = equationsToMatrix(str2sym(operation_info(1:3)), vars);

            rotation_mat = double(R_abc2xyz * rotation_mat * R_abc2xyz^(-1));
            
            if length(operation_info)==3
                symop_xyz = [symop_xyz, Oper(rotation_mat)];
            else
                is_conjugate = str2double(operation_info{4})==-1;
                symop_xyz = [symop_xyz, Oper(rotation_mat, 'conjugate', is_conjugate)];
            end
        end
    end
end
frewind(fid);

cif_struct.symop_xyz = symop_xyz(2:end);
%% magnetic cif
if ~contains(filename, "mcif")
    return
end

line_moment_loop_start = false;

while ~feof(fid)
    line = fgetl(fid);
    
    if contains(line, '_space_group_magn.transform_BNS_Pp_abc')
        line_transform = line;
    end
    
    if contains(line, '_space_group_magn.number_BNS')
        line_number_BNS = line;
    end

    if contains(line, '_atom_site_moment.label')
        line_moment_loop_start = true;
    end

    if line_moment_loop_start
        if line(1)=="_"
            continue
        end

        line_moment = split(line);
        moment = str2double(str_simplify(line_moment(2:4)));
        if sum(abs(moment))~=0
            easy_axis_abc = moment;
            easy_axis_xyz = R_abc2xyz * easy_axis_abc;
            line_moment_loop_start = false;
        end
    end
end
frewind(fid);

line_number_BNS = split(line_number_BNS);
number_BNS = split( regexprep(line_number_BNS{2}, '^[''"]|[''"]$', '') , '.'); % remove quotes

if ~isempty(line_transform) % assume it is from findsym and must be eye(3)
    line_transform = split(line_transform);
    transform_BNS_Pp_abc = regexprep(line_transform{2}, '^[''"]|[''"]$', ''); % remove quotes
    transform_BNS_P = split(transform_BNS_Pp_abc,";");
    transform_BNS_P = split(transform_BNS_P{1},',');

    transform_BNS_P = regexprep(transform_BNS_P, '(\d)([a-z])', '$1*$2'); % 2a -> 2*a
    
    syms a b c real;
    vars = [a,b,c];
    
    transform_BNS_P_mat = equationsToMatrix(str2sym(transform_BNS_P), vars);
else
    transform_BNS_P_mat = sym(eye(3));
end

cif_struct.number_BNS = [str2double(number_BNS{1}), str2double(number_BNS{2})];
cif_struct.easy_axis_abc = easy_axis_abc / norm(easy_axis_abc);
cif_struct.easy_axis_xyz = easy_axis_xyz / norm(easy_axis_xyz);
cif_struct.transform_BNS_P_mat = transform_BNS_P_mat;

fclose(fid);
end

function str2 = str_simplify(str1)
str2 = regexprep(str1, '\(.*', '');
end