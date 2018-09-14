function [ f_inv ] = numeric_inverse_x_min_sinh(xmin,xmax,coef,y )
%NUMERIC_INVERSE_X_MIN_SINH Summary of this function goes here
%   Gives a value x that satisfies y = x - coef(1) * sinh(coef(2)*x) 
% If the function is invertible, this function will return the only value y
% that satisfies the previous property if the solution is sufficiently
% close to the mean value of xmin and xmax (it does not need to be in the interval).
% If it is not invertible, it will look for a solution in
% the interval [xmin xmax].
problem.objective = @(x) (x + coef(1)*sinh(coef(2)*x - y));
%if coef(1)*coef(2) >= 1
%    problem.x0 = (xmin + xmax) / 2;
%else
    problem.x0 = [min(xmin,xmax),max(xmin,xmax)];
%end
problem.solver = 'fzero';
problem.options = optimset(@fzero);

f_inv = fzero(problem);

end

