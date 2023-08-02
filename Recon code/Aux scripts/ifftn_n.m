function y = ifftn_n(varargin)

if nargin>1 % ifftn(x,sz)
    factor = sqrt(prod(varargin{2}));
else % ifftn(x)
    factor = sqrt(numel(varargin{1}));
end

y = ifftn(varargin{1})*factor;

end