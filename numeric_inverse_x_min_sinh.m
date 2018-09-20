function [ f_inv ] = numeric_inverse_x_min_sinh(xmin,xmax,coef,y )
%NUMERIC_INVERSE_X_MIN_SINH Summary of this function goes here
%   Gives a value x that satisfies y = x - coef(1) * sinh(coef(2)*x) in the
% interval xmin,xmax, if there exists one.

warning('The function numeric_inverse_x_min_sinh is used, and it shouldn''t');

f_inv = zeros(1,length(y));

for i = 1:length(y)
    problem.objective = @(x) (x + coef(1)*sinh(coef(2)*x - y(i)));
    problem.x0 = [min(xmin,xmax),max(xmin,xmax)];
    problem.solver = 'fzero';
    problem.options = optimset(@fzero);
    f_inv(i) = fzero(problem)
end

end

