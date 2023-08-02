function x = CEadj(y,reconopts)
%
% x = CE(y,reconopts)
% 
% Adjoint operator in "continuous" space
%
% Input: y, data in k-space (column vector or any shape)
%
% Output: x, data in x-space (column vector)
%
% if data is uniform (reconopts.Uniform), it assumes y to be fftshift'ed
%
% reconopts: recon information:
%   - Nd
%   - Uniform
%   - rawdata_dim
%   - pix
%   - df
%   - t
%   - FTst
%   - w

Nd = reconopts.Nd;
pix = reconopts.pix;
df = reconopts.df;
t = reconopts.t;
if reconopts.Uniform
    t = fftshift(t);
    pix = fftshift(pix);
end

y = reshape(y,reconopts.rawdata_dim);
x = zeros(Nd);

if reconopts.Uniform
    for pj = min(pix(:)):max(pix(:))
        mask = (pix == pj);
        if sum(mask(:))>0 % It only does it if that frequency exists
            pp = pj*df;
            x = x + mask.*ifftn(y.*exp(1i*2*pi*t*pp));
        end
    end

    x = x(:)*sqrt(prod(Nd)); % Normalized in x-y (not in f)

else
    y = y.*reconopts.w;
    % For the NUFFT it is worth parallelizing?
    for pj = min(pix(:)):max(pix(:))
        mask = (pix == pj);
        if sum(mask(:))>0 % It only does it if that frequency exists
            pp = pj*df;
            x = x + mask.*nufft_adj(y.*exp(1i*2*pi*t*pp),reconopts.FTst).';
        end
    end

    x = x(:)/sqrt(prod(Nd));

end
