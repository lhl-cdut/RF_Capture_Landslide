function [ind,MAic] = M_aic(x)

L = length(x);
MAic = zeros(L, 1);
ind = 0;
for k = 1:1:L
    MAic(k) = k*log10(var(x(1:k)))+(L-k-1)*log10(var(x(k+1:L)));
end

MAic(MAic == 0) = [];
ind = find(MAic == min(MAic(2:end)));