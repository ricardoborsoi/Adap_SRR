% computes O-MOMS function at specified value(s) and degree
% currently degrees of 0 - 5 is supported

% changelog:
% 1/15/2010 Fixed a bug in degree = 4 case
% 4/10/2010 Handle x = 1/2 at degree = 2 separately. the value should be 1/2

% Copyright (c) 2010, Meng Wang
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are 
% met:
% 
% Redistributions of source code must retain the above copyright notice, 
% this list of conditions and the following disclaimer.
% Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution. Neither the name
% of CUHK nor the names of its contributors may be used to endorse or 
% promote products derived from this software without specific prior 
% written permission.
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS 
% IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
% THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
% PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
% CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

function value = omoms_at(x, degree)
value = zeros(size(x));
if (nargin < 2)
    fprintf('Need to specify degree.\n');
    return;
end;

absx = abs(x);

switch (degree)
    case 0
        value(absx < 1/2) = 1;
    case 1
        logical = absx < 1;
        value(logical) = 1 - absx(logical);
    case 2
        logical = absx < 1/2;
        value(logical) = -absx(logical).^2 + 43/60;
        logical = (absx >= 1/2 & absx < 3/2);
        value(logical) = absx(logical).^2/2 - 3*absx(logical)/2 + 137/120;
        logical = (absx == 1/2); % we handle 1/2 differently. the value ought to be 1/2
        value(logical) = 60/120;
    case 3
        logical = absx < 1;
        value(logical) = absx(logical).^3/2 - absx(logical).^2 + absx(logical)/14 + 13/21;
        logical = absx >= 1 & absx < 2;
        value(logical) = -absx(logical).^3/6 + absx(logical).^2 - 85/42*absx(logical) + 58/42;
    case 4
        logical = absx < 1/2;
        value(logical) = absx(logical).^4/4 - 13*absx(logical).^2/24 - absx(logical)/1512 + 11383/20160;
        logical = absx >= 1/2 & absx < 1;
        value(logical) = -absx(logical).^4/6+5*absx(logical).^3/6-47*absx(logical).^2/36+131*absx(logical)/378+1693/3360;
        logical = absx >= 1 & absx < 3/2;
        value(logical) = -absx(logical).^4/6+5*absx(logical).^3/6-47*absx(logical).^2/36+1051*absx(logical)/3024+5069/10080;
        logical = absx >= 3/2 & absx < 2;
        value(logical) = absx(logical).^4/24 - 5*absx(logical).^3/12 + 227*absx(logical).^2/144 - 2021*absx(logical)/756 + 69101/40320;
        logical = absx >= 2 & absx < 5/2;
        value(logical) = absx(logical).^4/24 - 5*absx(logical).^3/12 + 227*absx(logical).^2/144 - 20213*absx(logical)/7560 + 69133/40320;
        logical = absx >= 5/2 & absx < 3;
        value(logical) = 1/15120*(3-absx(logical));
    case 5
        logical = absx < 1/2;
        value(logical) = -absx(logical).^5/12 + absx(logical).^4/4 - 5*absx(logical).^3/99 -325*absx(logical).^2/792 + 749/1440;
        logical = absx >= 1/2 & absx < 1;
        value(logical) = -absx(logical).^5/12 + absx(logical).^4/4 - 5*absx(logical).^3/99 -431*absx(logical).^2/1056 - 7*absx(logical)/3168 + 10997/21120;
        logical = absx >= 1 & absx < 3/2;
        value(logical) = absx(logical).^5/24 - 3*absx(logical).^4/8 + 505*absx(logical).^3/396 - 181*absx(logical).^2/96 + 2693*absx(logical)/3168 + 6757/21120;
        logical = absx >= 3/2 & absx < 2;
        value(logical) = absx(logical).^5/24 - 3*absx(logical).^4/8 + 505*absx(logical).^3/396 - 4981*absx(logical).^2/2640 + 1691*absx(logical)/1980 + 3347/10560;
        logical = absx >= 2 & absx < 5/2;
        value(logical) = -absx(logical).^5/120 + absx(logical).^4/8 - 299*absx(logical).^3/396 + 6059*absx(logical).^2/2640 - 6949*absx(logical)/1980 + 691/320;
        logical = absx >= 5/2 & absx < 3;
        value(logical) = -absx(logical).^5/120 + absx(logical).^4/8 - 299*absx(logical).^3/396 + 36361*absx(logical).^2/15840 - 5057*absx(logical)/1440+136993/63360;
        logical = absx >= 3 & absx < 7/2;
        value(logical) = 1/15840*(7/2 - absx(logical)).^2;
    otherwise
        fprintf('Unsupported degree %d.\n', degree);
end;