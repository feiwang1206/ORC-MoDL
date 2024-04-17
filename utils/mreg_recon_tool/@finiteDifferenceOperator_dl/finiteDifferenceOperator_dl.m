function  A = finiteDifferenceOperator_dl(direction)

% usage:
%    A = finiteDifferenceOperator(direction)
%
% direction = array dimension in which the operation is applied

s.adjoint = 0;
s.direction = direction;

A = class(s,'finiteDifferenceOperator_dl');