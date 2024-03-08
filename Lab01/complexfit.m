function [a,b] = complexfit(Gk, n, uf, fs)
% truncate the data
w = uf/(fs/2)*pi  
% for i=1:N
%     if w(i) > wmax
%         Gk(i:end) = [];
%         break
%     end
% end

% calculate the e^(jnw) matrix
order = 1:n;
newN = length(Gk);
e_b = zeros(newN, n);
e_a = zeros(newN, n);
for i=1:newN
    e_b(i,:) = exp(-j*order*w(i));
end
% generate Y
for i=1:newN
    e_a(i,:)=e_b(i,:)*(-Gk(i)) ;
end
Y = [e_a e_b];
% generate Yaug , Xaug
Yaug = [real(Y);imag(Y)];
Xaug = [real(Gk);imag(Gk)];
theta = pinv(Yaug)*Xaug;
a = theta(1:n);
b = theta(n+1:end);
end