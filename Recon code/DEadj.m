function x = DEadj(y,reconopts)
%
% x = DEadn(y,reconopts)
% 
% Adjoint operator in discrete space
% x: data in x-space
% y: data in k-space
% reconopts: recon information:
%   - Nd,Nf
%   - Uniform
%   - rawdata_dim
%   - linpix
%   - lintix
%   - idxtix
%   - FTst
%   - w

Nd = reconopts.Nd;
Nf = reconopts.Nf;
linpix = reconopts.linpix;
lintix = reconopts.lintix;
idxtix = reconopts.idxtix;

y = reshape(y,reconopts.rawdata_dim);
M_hat = zeros([Nf reconopts.rawdata_dim]);

if reconopts.Uniform

    % I1) Augment m_hat(kx) to M_hat(kx,ky,kf) = m_hat(kx,ky) \delta( kf-t(kx,ky) )
    M_hat(lintix) = y;
    % I2) Use inverse FT to obtain adjoint
    M = ifftn_n(M_hat)*sqrt(Nf); % Normalized only in x-y (not f)

else

    % It separates space and frequency and do first space (more
    % efficient) and then frequency. There could be a small problem
    % with gridding with many zeroes.
    % I1) Augment m_hat(kx) to M_hat(kx,ky,kf) = m_hat(kx,ky) \delta( kf-t(kx,ky) )
    M_hat(lintix) = y.*reconopts.w;
    % I2) Inverse FT M_hat(kx,ky,kf) to obtain M(x,y,f)
    M1 = zeros([Nf Nd]);
    for f = idxtix
        M1(f,:,:) = nufft_adj(M_hat(f,:).',reconopts.FTst).'; % It shouldn't use dcf
    end
    M = ifft(M1,[],1)/sqrt(prod(Nd))*Nf;

end

% I3) Sample M to obtain mp(x,y) = M( x,y p(x,y) )
x = reshape(M(linpix),Nd);
x = x(:);
