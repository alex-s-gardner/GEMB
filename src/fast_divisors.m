function d = fast_divisors(n)
% fast_divisors is a faster implementation of the divisors function and
% does NOT require the Symbolic Math Toolbox. 
% 
%% Syntax
% 
%  d = fast_divisors(n)
%
%% Description
% 
% d = fast_divisors(n) finds all nonnegative divisors of an integer n.
% 
%% Example 
% 
%  fast_divisors(42)
%  ans =
%     1     2     3     6     7    14    21    42
%
%% Author Information
% The Glacier Energy and Mass Balance (GEMB) was created by Alex Gardner, with contributions
% from Nicole-Jeanne Schlegel and Chad Greene. Complete code and documentation are available
% at https://github.com/alex-s-gardner/GEMB. Please cite any use of GEMB as:
% 
% Gardner, A. S., Schlegel, N.-J., and Larour, E.: Glacier Energy and Mass Balance (GEMB): 
% a model of firn processes for cryosphere research, Geosci. Model Dev., 16, 2277–2302, 
% https://doi.org/10.5194/gmd-16-2277-2023, 2023. 

% Find divisors up to the square root
k = 1:ceil(sqrt(n));
d = k(rem(n,k)==0);

% Find corresponding divisors > sqrt(N) and combine
d = unique([d n./d]);

end