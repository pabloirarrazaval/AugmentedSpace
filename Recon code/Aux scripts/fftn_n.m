function y = fftn_n(varargin)

if nargin>1 % fftn(x,sz)
    factor = sqrt(prod(varargin{2}));
else % fftn(x)
    factor = sqrt(numel(varargin{1}));
end

y = fftn(varargin{1})/factor;

end