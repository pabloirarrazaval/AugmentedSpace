function st = prepares_nufft(data)

% Check that previous functions were called
if isempty(data.kx) || isempty(data.ky)
    error('Trajectory not defined (fill in data.kx and data.ky).'); end

% Prepare structure to call NUFFT
Nd = data.Nd;
om = [data.kx data.ky]*2*pi; % From -pi to pi
Jd = [5 5];       % Kernel extension
Kd = Nd*2;        % Oversampling (times 2)
n_shift = data.shift;   % Origin at center when data.shift = data.Nd/2

st = nufft_init(om,Nd,Jd,Kd,n_shift);

end