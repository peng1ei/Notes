function [r] = dirichlet_rnd(a,dim)

%DIRICHLET_RND Random matrices from dirichlet distribution.
%   R = DIRICHLET_RND(A,DIM) returns a matrix of random numbers chosen   
%   from the dirichlet distribution with parameters vector A.
%   Size of R is (N x N) where N is the size of A or (N x DIM) if DIM is given.

if nargin < 1, error('Requires at least one input arguments.'); end
if nargin > 2, error('Requires at most two input arguments.'); end

[rows columns] = size(a); 
if nargin == 1, dim = rows * columns; end
if nargin == 2, 
   if prod(size(dim)) ~= 1, error('The second parameter must be a scalar.'); end
end
if rows~=1 & columns~=1, error('Requires a vector as an argument.'); end
if any( a < 0 | ~isreal(a)), error('Parameters of Dirichlet function must be real positive values.');end

% fastest method that generalize method used to sample from
% the BETA distribuition: Draw x1,...,xk from independent gamma 
% distribuitions with common scale and shape parameters a1,...,ak, 
% and for each j let rj=xj/(x1+...+xk).

N = rows * columns;
for i = 1 : N
    x(:,i) = gamrnd(a(i),a(i),dim,1); % generates dim random variables 
end                                   % with gamma distribution
r = x./repmat(sum(x,2),[1 N]);

% For more details see "Bayesian Data Analysis" Appendix A




