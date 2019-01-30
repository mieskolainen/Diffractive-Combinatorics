% Poincare disk distance
%
% input  u = n-vector (l2-norm < 1)
%        V = (M x n) matrix of n-vectors (l2-norm < 1)
%
% mikael.mieskolainen@cern.ch, 2018
% Licensed under the MIT License <http://opensource.org/licenses/MIT>.

function D = poincaredist(u,V)

D = zeros(size(V,1),1);
%
% Poincare disk distance
for i = 1:size(D,1)
    v = V(i,:);   
    D(i) = acosh(1 + 2*norm(u-v)^2/((1-norm(u)^2)*(1-norm(v)^2)) );
end
%}
%{
% Euclidian
for i = 1:size(D,1)
    v = V(i,:);
    D(i) = norm(u - v);
end
%}

end