function d = pixtix(d,op)
%
% data = pixtix(data,operation)
%
% Create or changes the values of pix and tix. These are the field map and 
% time map as indexes to the matrices (changes Hz or s to pixels).
%
% operation = 'create +':   Creates them in natural order in all dimensions
%           = 'create -':   Creates them in natural order (simple
%                           discretization in f, includes negatives)
%           = 'shift x':    fftshift's in the image/k-space domain
%           = 'shift f':    fftshift's in the frequency/tme domain
%           = 'linear':     Creates linpix and lintix. Linearized indices
%                           to arrays.
%           = 'avail':      Creates idxpix and idxtix. Indices of pix and
%                           tix with data. To avoid computing fft's for
%                           only zeros.

Nd = d.Nd;

switch op
    case 'create +'
        d.pix = round(d.p/d.df) + d.Nf/2 + 1;
        d.tix = round(d.t/d.dt) + d.Nf/2 + 1;

    case 'create -'
        d.pix = round(d.p/d.df);
        d.tix = round(d.t/d.dt);

    case 'shift x'
        d.pix = fftshift(d.pix);
        d.tix = fftshift(d.tix);

    case 'shift f'
        d.pix = d.pix - d.Nf/2 + (d.pix < (d.Nf/2+1)).*d.Nf;
        d.tix = d.tix - d.Nf/2 + (d.tix < (d.Nf/2+1)).*d.Nf;

    case 'linear'
        switch length(Nd)
            case 1 % Only cartesians
                % Transforms p and t (kf) to indexes into the range -Nf/2 to Nf/2-1
                d.linpix = ((0:Nd-1).')*d.Nf + d.pix(:);
                d.lintix = ((0:Nd-1).')*d.Nf + d.tix(:);
            case 2
                [xx,yy] = meshgrid(0:Nd(2)-1,0:Nd(1)-1);
                d.linpix = xx(:)*Nd(1)*d.Nf + yy(:)*d.Nf + d.pix(:);
                if d.Uniform
                    d.lintix = xx(:)*Nd(1)*d.Nf + yy(:)*d.Nf + d.tix(:); 
                else
                    d.lintix = (0:d.L-1).'*d.Nf + d.tix(:);
                end
        end

    case 'avail'
        availpix = zeros([1 d.Nf]); availtix = zeros([1 d.Nf]);
        for ii = 1:d.Nf
            availpix(ii) = (sum(d.pix(:)==ii)>0);
            availtix(ii) = (sum(d.tix(:)==ii)>0);
        end
        d.idxpix = find(availpix>0);
        d.idxtix = find(availtix>0);
        
    otherwise
        error('Operation %s not defined.',op);
        
end