function ndx = my_sub2ind(siz,sub)

siz = double(siz);
k = [1 cumprod(siz(1:end-1))];
ndx = 1;
for i = 1:numel(siz)
    v = sub(:,i);
    ndx = ndx + (v-1)*k(i);
end
