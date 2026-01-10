function D = fast_divisors(N)
% fast_divisors is a faster implementation of the divisors function and
% does NOT require the Symbolic Math Toolbox. 
% 
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
K = 1:ceil(sqrt(N));
D = K(rem(N,K)==0);

% Find corresponding divisors > sqrt(N) and combine
D = [D sort(N./D)];

end