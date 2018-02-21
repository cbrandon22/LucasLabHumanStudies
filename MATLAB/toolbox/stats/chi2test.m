function[p,chi2,df] = chi2test(counts)

columnSums = sum(counts,1);
rowSums = sum(counts,2);

sumSum = sum(columnSums);

expected = zeros(size(counts));
for i = 1:size(counts,1)
    for j = 1:size(counts,2)
        expected(i,j) = rowSums(i)*columnSums(j)/sumSum;
    end
end

df = (size(counts,1)-1)*(size(counts,2)-1);

chi2 = sum(sum(((counts - expected).^2)./expected));

p = 1-chi2cdf(chi2,df);