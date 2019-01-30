% Binary to Gray code
%
% mikael.mieskolainen@cern.ch, 2018
% Licensed under the MIT License <http://opensource.org/licenses/MIT>.

function g = bin2gray(b)
g(1) = b(1);
for i = 2:length(b)
    x = xor(str2num(b(i-1)), str2num(b(i)));
    g(i) = num2str(x);
end

end