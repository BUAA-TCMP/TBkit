function TBkitobj = input_orb_struct(TBkitobj,filename,mode,options)
%INPUT_ORB_STRUCT Import crystal structure data and populate TBkit object
%   This function reads structural information from POSCAR-like files and 
%   populates the TBkit object with orbital, symmetry, and quantum number data
%
%   Inputs:
%       TBkitobj  - Target TBkit object to populate
%       filename  - Input file name (default: 'POSCAR')
%       mode      - Input format mode: 'vasp', 'tbsk', or 'sym' (default: 'vasp')
%       options   - Structure containing processing options:
%           symbolic    : Enable symbolic precision (default: false)
%           Operation   : Enable symmetry operations (default: false)
%           warning     : Show warnings (default: true)
%           spin        : Spin treatment mode: 'spinless', 'wannier', 'block'
%
%   Output:
%       TBkitobj - Updated TBkit object with structural data
%
%   Example:
%       TBobj = input_orb_struct(TBobj, 'myPOSCAR', 'tbsk', ...
%                struct('symbolic',true,'spin','wannier'))

arguments
    TBkitobj TBkit;
    filename string = 'POSCAR';
    mode char {mustBeMember(mode,{'vasp','tbsk','sym','s,sz'})} = 'vasp';
    options.symbolic logical = false;
    options.Operation logical = false;
    options.warning logical = true;
    options.spin char {mustBeMember(options.spin,{'spinless','wannier','block'})} = 'spinless';
end

% Initialize symmetry library if needed
if options.Operation
    import spglib_matlab.*;
end

% Read basic structure information from file
[TBkitobj.Rm, tmpsites, TBkitobj.Atom_name, TBkitobj.Atom_num, elements] = POSCAR_read(filename, mode);

% Validate atomic site counts
if options.warning
    validate_site_counts(TBkitobj, tmpsites);
end

% Process sites based on input mode
switch mode
    case 'vasp'
        % Direct assignment for VASP format
        TBkitobj.sites = tmpsites;
        
    case {'tbsk','sym'}
        % Generate complete site information with spin consideration
        tmpsites = TBkit.MakeUpSites(tmpsites, options.spin);
end

% Store orbital positions in lattice coordinates
TBkitobj.orbL = [[tmpsites.rc1].', [tmpsites.rc2].', [tmpsites.rc3].'];
TBkitobj.orbL_cart = TBkitobj.orbL * TBkitobj.Rm;

% Process symmetry operations if enabled
if options.Operation
    TBkitobj = process_symmetry_operations(TBkitobj);
end

% Process quantum numbers and orbital symmetries
if ~isempty(elements)
    elements.Properties.RowNames = elements.atom_symbol;
end

switch mode
    case 'tbsk'
        TBkitobj = process_tbsk_mode(TBkitobj, tmpsites, elements);
    case {'sym','s,sz'}
        TBkitobj = process_sym_mode(TBkitobj, tmpsites, elements);
end

% Handle symbolic precision conversion
if options.symbolic
    TBkitobj = handle_symbolic_conversion(TBkitobj);
end

end
%% Helper functions
function validate_site_counts(TBkitobj, tmpsites)
% VALIDATE_SITE_COUNTS Check consistency of atomic site counts
if sum(TBkitobj.Atom_num) ~= TBkitobj.Basis_num
    if isa(TBkitobj,'HR')
        objType = TBkitobj.Basis_num;
        warning(['Sum of POSCAR 7th line differs from %s.\n'...
            'Using %s enforcedly'], sum(TBkitobj.Atom_num), objType);

    else
    end
end

if length(tmpsites) ~= TBkitobj.Basis_num
    if isa(TBkitobj,'HR')
        objType = TBkitobj.Basis_num;
        warning(['Sum of POSCAR 7th line differs from %s.\n'...
            'Using %s enforcedly'], length(tmpsites), objType);

    else
    end
end
end

function TBkitobj = process_symmetry_operations(TBkitobj)
% PROCESS_SYMMETRY_OPERATIONS Calculate space group symmetry operations
types = arrayfun(@(i) ones(1,TBkitobj.Atom_num(i))*i, ...
    1:length(TBkitobj.Atom_num), 'UniformOutput', false);
types = [types{:}];

spglib_database = spglib_matlab.spg_get_dataset_from_sites(...
    TBkitobj.Rm, TBkitobj.orbL.', types);

TBkitobj.sgn = spglib_database.spacegroup_number;
TBkitobj.SymmetryOperation = [spglib_database.rotations; 
                             spglib_database.translations];
end

function TBkitobj = process_tbsk_mode(TBkitobj, tmpsites, elements)
% PROCESS_TBSK_MODE Process orbital data for TBSK format
tmp_orb_symL = [tmpsites.orb_sym].';
TBkitobj.Basis_num = size(TBkitobj.orbL, 1);

for i = 1:TBkitobj.Basis_num
    element_data = table2array(elements(...
        char(tmpsites(i).element), {'atom_number','n'}));
    
    TBkitobj.elementL(i,:) = element_data(1);
    TBkitobj.quantumL(i,:) = process_quantum_numbers(...
        element_data(2), tmpsites(i), tmp_orb_symL(i,:));
    
    try
        TBkitobj.orb_symL(i,:) = TBkit.Ymlsym(...
            TBkitobj.quantumL(i,2), TBkitobj.quantumL(i,3), ...
            tmp_orb_symL(i,:));
    catch
        % Handle orbital symmetry conversion errors
    end
end
TBkitobj.sites = tmpsites;
end

function TBkitobj = process_sym_mode(TBkitobj, tmpsites, elements)
% PROCESS_SYM_MODE Process symmetry-related quantum numbers
tmp_orb_symL = [tmpsites.orb_sym].';

for i = 1:TBkitobj.Basis_num
    element_data = table2array(elements(...
        char(tmpsites(i).element), {'atom_number','n'}));
    
    TBkitobj.elementL(i,:) = element_data(1);
    quantum_numbers = process_quantum_numbers(...
        element_data(2), tmpsites(i), tmp_orb_symL(i,:));
    
    % Store additional J quantum number for symmetry mode
    TBkitobj.quantumL(i,1:4) = quantum_numbers;
    try
        TBkitobj.quantumL(i,5) = tmpsites(i).Jz;
    catch
    end
    
    TBkitobj.orb_symL(i,:) = TBkit.Ymlsym(...
        quantum_numbers(2), quantum_numbers(3), tmp_orb_symL(i,:));
end
end

function qnums = process_quantum_numbers(n, site, orb_symbk)
% PROCESS_QUANTUM_NUMBERS Calculate quantum numbers from site data
l = TBkit.orb2l(site.orb);
m = TBkit.orb_sym2m(orb_symbk);
if isfield(site,'spin')
    qnums = [n, l, m, site.spin];
else
    qnums = [n, l, m, nan];
end
end

function TBkitobj = handle_symbolic_conversion(TBkitobj)
% HANDLE_SYMBOLIC_CONVERSION Convert numerical values to symbolic format
precision_warning = false;

% Process lattice vectors
TBkitobj.Rm = check_symbolic_precision(TBkitobj.Rm, 'Rm', 1e-13);
precision_warning = precision_warning || any(cellfun(@length, ...
    arrayfun(@(x) char(string(x)), TBkitobj.Rm(:), 'UniformOutput', false)) > 15);

% Process orbital positions
TBkitobj.orbL = check_symbolic_precision(TBkitobj.orbL, 'orbL', 1e-13);
precision_warning = precision_warning || any(cellfun(@length, ...
    arrayfun(@(x) char(string(x)), TBkitobj.orbL(:), 'UniformOutput', false)) > 15);

if precision_warning
    display_precision_warning();
end
end

function sym_array = check_symbolic_precision(num_array, array_name, tol)
% CHECK_SYMBOLIC_PRECISION Validate numerical precision for conversion
sym_array = sym(num_array);
for idx = 1:numel(sym_array)
    num_val = num_array(idx);
    sym_val = sym_array(idx);
    
    if abs(double(sym_val) - num_val) > tol
        [i,j] = ind2sub(size(num_array), idx);
        warning('Precision loss in %s(%d,%d): %s', ...
            array_name, i, j, char(sym_val));
    end
end
end

function display_precision_warning()
% DISPLAY_PRECISION_WARNING Show precision recommendation
fprintf(['\n************* WARNING *****************\n'...
    'Numerical precision should reach 1e-13. Examples:\n'...
    '0.333333333333333 -> %s\n'...
    '0.866025403784439 -> %s\n'...
    'Consider using ''format long'' to verify input precision\n'],...
    string(sym(0.333333333333333)), string(sym(0.866025403784439)));
end