function y = fft_n(varargin)

% find first dimension that is not of size one
sz = size(varargin{1});
ix = 1;
while ix <= length(sz) && sz(ix) == 1
    ix = ix+1;
end
factor = sqrt(size(varargin{1},ix));

if nargin>1 && ~isempty(varargin{2}) % fft(x,n,dim|[])
    factor = sqrt(varargin{2});
elseif nargin>2 % fft(x,[],dim)
    factor = sqrt(size(varargin{1},varargin{3}));
end

y = fft(varargin{1})/factor;

end