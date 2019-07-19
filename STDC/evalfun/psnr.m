function result = psnr(X,Y)
X = double(X); 
Y = double(Y); 
m = mean((X(:)-Y(:)).^2); 
peak = max(X(:));
result = 10*log10(peak^2/m);