function [ f, fErr ] = propError(f, varlist, vals, errs )
% [f, fErr] = propError('f', {'var1', 'var2', ...}, [val1, val2, ...], {val1Err, val2Err, ...})
% calcualtes the function f and its error according to partial
% derivatives.
% f is a string, e.g. 'a*b';
% varlist is a cell array of strings, e.g {'a', 'b'}
% vals is a cell array of vectors of values, e,g, {[1 2 3], [4 5 6]}
% errs is a cell array of vectors of error per parameter, e.g. {[0.1 0.21 0.05], [0.07 0.06 0.01]}

% Yair Margalit 12/11/2016
%original function:
%(c) Brad Ridder 2007. Feel free to use this under the BSD guidelines. If
%you wish to add to this program, just leave my name and add yours to it.

vals = cellfun(@makeColumn, vals, 'UniformOutput', 0); %make sure input is column vectors
errs = cellfun(@makeColumn, errs, 'UniformOutput', 0); %make sure input is column vectors
n = numel(varlist);
for i = 1 : n %transform strings to syms
    syms(varlist{i})
end
f = eval(f); %transform f from string to symbolic
sig = vpa(ones(1,n));
for i = 1:n %calcualte partial derivatives
    sig(i) = diff(f,varlist{i},1);
end
% errs = repmat(errs, length(vals{1}), 1);
fErr = sqrt( sum( ( subs(sig,varlist,vals) .* cell2mat(errs) ).^2, 2) ); %error in f according to partial derivatives
fErr = double(fErr);
f = double(subs(f,varlist,vals));

% sigma = [{subs(f,varlist,vals)} {'+/-'} {error};
%          {'Percent Error'} {'+/-'} {abs(100*(error)/subs(f,varlist,vals))}];

end

function [v] = makeColumn(v)
if isrow(v)
    v=v';
end
end

