function C = lt(A,B)
%LT Less than operation for HR objects
%   Specialized comparison/operation function for HR objects.
%
%   Syntax:
%       C = lt(A,B)
%
%   Inputs:
%       A - HR object
%       B - HR object or other compatible type
%
%   Outputs:
%       C - Result of operation
%
%   Note:
%       Currently supports limited operations between HR objects and
%       other data types.
%
%   Example:
%       result = lt(H_hr, 'POSCAR'); % Load POSCAR into HR object
if isa(A,'HR') && isa(B,'HR')
H_hr1 = A;
H_hr2 = B;
if H_hr1.WAN_NUM ~= H_hr2.WAN_NUM
error('WAN_NUM different');
end
error('not support at present.')
elseif isa(A,'HR') && ~isa(B,'HR')
switch class(B)
case 'string'
C = lt(A,char(B));
return;
case 'char'
switch B(1)
case {'P','p'}
C = A.input_orb_struct(B,'tbsk');
case {'w','W'}
case {'k','K'}
C = A.kpathgen3D(B);
otherwise
end
case 'double'
switch size(B,1)
case A.WAN_NUN
C = A.input_orb_init(B);
otherwise
end
otherwise
end
elseif ~isa(A,'HR') && isa(B,'HR')
error('not support at present.');
end
end
