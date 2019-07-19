function sub = my_ind2sub(siz,ndx)

siz = double(siz);
sub = zeros(numel(siz),numel(ndx));
k = [1 cumprod(siz(1:end-1))];
for i = numel(siz):-1:1,
    vi = rem(ndx-1, k(i)) + 1;
    vj = (ndx - vi)/k(i) + 1;
    sub(i,:) = vj;
    ndx = vi;
end
