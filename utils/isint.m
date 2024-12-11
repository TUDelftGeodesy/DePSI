function res = isint(in)
% ISINT  True for integers.
%   ISINT(IN) returns 1 for those elements that are integers,
%   and 0 otherwise.  Valid for numbers in [-2147483648, 2147483647].
%   Mainly used for input checking in functions.
%
%   To test if any/all element of an array is an integer, use any or all:
%     if(any(isint(1:4))) disp('at least one element is an integer'); end;
%     if(all(isint(1:4))) disp('all elements are integers'); end;
%
%   Example:
%     a = [1 2 3 4 5.5 8.34];
%     isint(a)%  returns [1 1 1 1 0 0]
%
%   See also INT32
%

% $Revision: 1.3 $  $Date: 2001/09/28 14:24:32 $
% Bert Kampes, 15/03/00

%res = 0;
%if (int32(in)==in) res = 1; end;

%// BK 07-Aug-2001
%%%nico sneeuw: res=(rem(in,1));
res=(int32(in)==in);

%%% EOF
