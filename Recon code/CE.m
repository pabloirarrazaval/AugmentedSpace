function y = CE(x,reconopts)
%
% y = CE(x,reconopts)
% 
% Acquisition operator in "continuous" space
%
% Input: x, data in x-space (column vector or any shape)
%
% Output: y, data in k-space (column vector)
%
% if data is uniform (reconopts.Uniform), it assumes x to be fftshift'ed
%
% reconopts: recon information:
%   - Nd
%   - Uniform
%   - rawdata_dim
%   - tix
%   - dt
%   - p
%   - FTst
%   - w

Nd = reconopts.Nd;
tix = reconopts.tix;
dt = reconopts.dt;
p = reconopts.p;
if reconopts.Uniform
    p = fftshift(p);
    tix = fftshift(tix);
end

x = reshape(x,Nd);
y = zeros(reconopts.rawdata_dim);

if reconopts.Uniform
    for tj = min(tix(:)):max(tix(:))
        mask = (tix == tj);
        if sum(mask(:))>0 % It only does it if that time exists
            tt = tj*dt;
            y = y + mask.*fftn(x.*exp(-1i*2*pi*tt*p));
        end
    end
    y = y(:)/sqrt(prod(Nd)); % Normalized in x-y (not in f)

else
    % For the NUFFT is it worth parallelizing?
    for tj = min(tix(:)):max(tix(:))
        mask = (tix == tj);
        if sum(mask(:))>0 % It only does it if that time exists
            tt = tj*dt;
            % nufft puts x first, needs to transpose
            y = y + mask.*nufft((x.*exp(-1i*2*pi*tt*p)).',reconopts.FTst);
        end
    end
    y = y(:).*reconopts.w(:)/sqrt(prod(Nd));

end
