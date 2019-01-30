% Create binary combinatorics matrix
%
% Default bit order (0,0), (1,0), (0,1), (1,1)
%
% input:      N  =  Number of elements
%       reverse  =  Reverse binary order (default false)
% output:   CBM  =  Binary matrix
%
% mikael.mieskolainen@cern.ch, 2018
% Licensed under the MIT License <http://opensource.org/licenses/MIT>.

function CBM = createCBM(N, reverse)

if (nargin == 1)
    reverse = false;
end

CBM = zeros(2^N, N);

% Create binary codes
for i = 0:size(CBM,1)-1
    binary_str = dec2bin(i, N);
    binary_vect = zeros(1, N);
    
    % dow mark ones
    for j = N:-1:1
       binary_vect(j) = str2double(binary_str(j)); 
    end
    CBM(i+1,:) = binary_vect;
end

CBM = fliplr(CBM);

% Reverse bit direction
if (reverse)
   CBM = fliplr(CBM); 
end

end