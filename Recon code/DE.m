function y = DE(x,reconopts)
%
% y = DE(x,reconopts)
% 
% Acquisition operator in discrete space
% x: data in x-space
% y: data in k-space
% reconopts: recon information:
%   - Nd,Nf,L
%   - Uniform
%   - rawdata_dim
%   - linpix
%   - lintix
%   - idxpix
%   - FTst
%   - w

Nd = reconopts.Nd;
Nf = reconopts.Nf;
linpix = reconopts.linpix;
idxpix = reconopts.idxpix;
lintix = reconopts.lintix;
L = reconopts.L;

x = reshape(x,Nd);

% D1) Augment m(x,y) to M(x,y,f) = m(x,y) \delta( whof-p(x,y) )
M = zeros([Nf Nd]);
M(linpix) = x;

% D2) FT M(x,y,f) to obtain M_hat(kx,ky,kf)
if reconopts.Uniform
    M_hat = fftn_n(M)*sqrt(Nf); % Normalized only in x-y (not f)
    % D3) Sample M_hat to obtain m_hat(k) = M_hat(k, t(k) )
    y = M_hat(lintix);
else
    % It separates space and frequency and do first space (more
    % efficient) and then frequency. There could be a small problem
    % with gridding with many zeroes.
    M1 = zeros([Nf L]);
    for f = idxpix
        mx = reshape(M(f,:),Nd).'; % nufft puts x first, needs to transpose
        M1(f,:) = nufft(mx,reconopts.FTst);
    end
    M_hat = fft(M1,[],1)/sqrt(prod(Nd));
    % D3) Sample M_hat to obtain m_hat(k) = M_hat(k, t(k) )
    y = M_hat(lintix).*reconopts.w;
end

