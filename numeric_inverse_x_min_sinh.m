function [ f_inv ] = numeric_inverse_x_min_sinh(xmin,xmax,coef,y )
%NUMERIC_INVERSE_X_MIN_SINH Summary of this function goes here
%   Gives a value x that satisfies y = x - coef(1) * sinh(coef(2)*x) in the
% interval xmin,xmax, if there exists one.

problem.objective = @(x) (x + coef(1)*sinh(coef(2)*x - y));
problem.x0 = [min(xmin,xmax),max(xmin,xmax)];
problem.solver = 'fzero';
problem.options = optimset(@fzero);

f_inv = fzero(problem);

end

