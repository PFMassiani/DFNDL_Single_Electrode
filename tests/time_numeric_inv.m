tic;
N = 10;
coefs = [1 2];
xmin = 0;
xmax = 7;
for i = 1:N
    numeric_inverse_x_min_sinh(i/1.5,xmin,xmax,coefs);
   p = i;
end
toc;