function varargout = TBSK_Var_gen_single(L_1,L_2,m_1,m_2,nn_level,l,m,n,options)
%TBSK_VAR_GEN_SINGLE Generate Slater-Koster variables for single orbital pair
%
%   [COEFF, V_STR, S_STR] = TBSK_VAR_GEN_SINGLE(L_1, L_2, M_1, M_2, NN_LEVEL, L, M, N)
%   Generates Slater-Koster coefficients and variable names for a single orbital pair
%
%   Inputs:
%       L_1, L_2   - Angular momentum quantum numbers (0-3)
%       m_1, m_2   - Magnetic quantum numbers (-L to L)
%       nn_level   - Nearest neighbor level (default = -1)
%       l, m, n    - Direction cosines (default = 0)
%       options.overlap - Include overlap terms (default = false)
%
%   Outputs:
%       Coeff      - Slater-Koster coefficients [3×1]
%       V_str      - Hopping variable names {'VssS', 'VssP', 'VssD'}
%       S_str      - Overlap variable names (when options.overlap=true)
%
%   Note:
%       Automatically handles L_1 > L_2 case with proper sign changes
%
%   See also TBSK_VAR_GEN, HR.TBSK_Coeff_gen
arguments
L_1 double{mustBeInteger};
L_2 double{mustBeInteger};
m_1 double{mustBeInteger};
m_2 double{mustBeInteger};
nn_level double{mustBeInteger} = -1;
l  = 0;
m  = 0;
n  = 0;
options.overlap logical = false;
end
if nargin < 6 || (isa(l,'sym') || isa(m,'sym') ||isa(n,'sym')  )
sym_mode = true;
if ~isa(l,'sym')
syms l real;
end
if ~isa(m,'sym')
syms m real;
end
if ~isa(n,'sym')
syms n real;
end
else
sym_mode = false;
end
if L_1 > L_2
if nargout == 1
Coeff = HR.TBSK_Var_gen_single(L_2,L_1,m_2,m_1,nn_level,-l,-m,-n,'overlap',options.overlap);
varargout{1} = Coeff;
return;
elseif nargout == 2
[E_hop,E_int] = ...
HR.TBSK_Var_gen_single(L_2,L_1,m_2,m_1,nn_level,-l,-m,-n,'overlap',options.overlap);
varargout{1} = E_hop;
varargout{2} = E_int;
return;
elseif nargout == 3
[Coeff,varargout{2},varargout{3}] = ...
HR.TBSK_Var_gen_single(L_2,L_1,m_2,m_1,nn_level,-l,-m,-n,'overlap',options.overlap);
varargout{1} = Coeff;
return;
end
end
if sym_mode
Coeff = HR.TBSK_Coeff_gen(L_1,L_2,m_1,m_2,l,m,n,'sym_mode',true);
else
Coeff = HR.TBSK_Coeff_gen(L_1,L_2,m_1,m_2,l,m,n);
end
if nargout == 1 &&  options.overlap
varargout{1} = Coeff;
return
else
switch L_1
case 0
Char1 = 's';
case 1
Char1 = 'p';
case 2
Char1 = 'd';
case 3
Char1 = 'f';
end
switch L_2
case 0
Char2 = 's';
case 1
Char2 = 'p';
case 2
Char2 = 'd';
case 3
Char2 = 'f';
end
E_hop{1} = ['V',Char1,Char2,'S'] ;
E_hop{2} = ['V',Char1,Char2,'P'] ;
E_hop{3} = ['V',Char1,Char2,'D'] ;
if options.overlap || nargout ~= 1
E_integral{1} = ['S',Char1,Char2,'S'] ;
E_integral{2} = ['S',Char1,Char2,'P'] ;
E_integral{3} = ['S',Char1,Char2,'D'] ;
else
E_integral = [];
end
if nn_level > 0
E_hop{1} = [E_hop{1},'_',num2str(nn_level)];
E_hop{2} = [E_hop{2},'_',num2str(nn_level)];
E_hop{3} = [E_hop{3},'_',num2str(nn_level)];
if options.overlap || nargout ~= 1
E_integral{1} = [E_integral{1},'_',num2str(nn_level)];
E_integral{2} = [E_integral{2},'_',num2str(nn_level)];
E_integral{3} = [E_integral{3},'_',num2str(nn_level)];
end
end
end
if nargout == 3
varargout{1} = Coeff;
varargout{2} = E_hop;
varargout{3} = E_integral;
return
elseif nargout == 2
E_hop_sym = [sym(E_hop{1},'real'),...
sym(E_hop{2},'real'),...
sym(E_hop{3},'real') ]*Coeff.';
E_integral_sym = [sym(E_integral{1},'real'),...
sym(E_integral{2},'real'),...
sym(E_integral{3},'real') ]*Coeff.';
varargout{1} = E_hop_sym;
varargout{2} = E_integral_sym;
elseif nargout == 1 && ~options.overlap
E_hop_sym = [sym(E_hop{1},'real'),...
sym(E_hop{2},'real'),...
sym(E_hop{3},'real') ]*Coeff.';
varargout{1} = E_hop_sym;
end
end
