function m = FSEGrecon(m_hat,data,varargin)
% Frequency segmented reconstruction.
%
% Cartesian or non-uniform sampling
%
% m = FSEGrecon(m_hat,data,[Interp])
%       m_hat:    acquired data (hat because it is in k_space)
%       data:     information about the problem
%                   .Nd: dimensions
%                   .Nf: maximum number of frequency bins
%                   .targetdf: desired resolution in frequency
%                   .t:  time map in seconds as a function of kx in fft order
%                   .FTst: NUFFT structure used in Non-uniform
%       Interp:   interpolation: 'none' (default), 'linear', 'MFI', 'fix'
%
% m = FSEGrecon(m_hat,data,'fix',[f0])
%                   When Interp = 'fix', it will demodulate at a fixed frequency
%                   given by f0 (defaults to 0 if not specified).
%
% m = FSEGrecon(m_hat,data,'MFI',[filename])
%                   When Interp = 'MFI', the optional 'filename', if exists,
%                   will be read with precomputed weights. If does not
%                   exist, it will be created.
% 
% it will demodulate at a fixed frequency
%                   given by f0 (defaults to 0 if not specified).

global VERBOSE

% Starts by getting number of dimensions and size
Nd = data.Nd;
DIM = length(Nd);
if DIM>2, error('It only works for 1D or 2D.'); end
if DIM==1
    Nd = [Nd 1]; % fix dimensions for 1D
    t = data.t(:); % makes sure t is along first dimension
else
    t = data.t;
end

% Check for optional parameters
Interpolation = 'none';
if nargin>2
    Interpolation = varargin{1};
    if strcmp(Interpolation,'fix')
        Ffix = 0;
        if nargin>3
            Ffix = varargin{2};
            if ~isnumeric(Ffix)
                error('f0 must be a number.');
            end
        end
    elseif strcmp(Interpolation,'MFI')
        filename = 'MFIPrecomputedWeights_generic';
        if nargin>3
            filename = varargin{2};
            if ~(isa(filename,'char')||isa(filename,'string'))
                error('filename must be a string.');
            end
        end
    end
end      

switch Interpolation

    %----------------------------------------------------------
    case 'fix' % Demodulates to one frequency, given by Ns

        if data.Uniform
            m = fftshift(ifftn_n(fftshift((m_hat.*exp(1i*2*pi*t*Ffix)))));
        else
            m = nufft_adj(m_hat.*exp(1i*2*pi*t(:)*Ffix).*data.dcf,data.FTst).'/sqrt(prod(data.Nd));
        end

        if VERBOSE, fprintf('FSEG: One fixed demodulation at %f\n',Ffix);
        end
        
    %----------------------------------------------------------
    case 'none' % Traditional Frequency Segmented

        % Converts field to indices into the matrices
        data.pix = round(data.p/data.df);

        ncount = 0;
        m = zeros(Nd);
        if data.Uniform
            for p_j = min(data.pix(:)):max(data.pix(:))
                mask = (data.pix == p_j);
                if sum(mask(:))>0 % It only does it if that frequency exists
                    pp = p_j*data.df;  % will demodulate to p(p_j)
                    m_j = fftshift(ifftn_n(fftshift((m_hat.*exp(1i*2*pi*t*pp)))));
                    m = m + mask.*m_j;
                    ncount = ncount+1;
                end
            end
        else
            for p_j = min(data.pix(:)):max(data.pix(:))
                mask = (data.pix == p_j);
                if sum(mask(:))>0 % It only does it if that frequency exists
                    pp = p_j*data.df; % will demodulate to p(p_j)
                    m_j = nufft_adj(m_hat.*exp(1i*2*pi*t(:)*pp).*data.dcf,data.FTst).';
                    m = m + mask.*m_j;
                    ncount = ncount+1;
                end
            end
            m = m/sqrt(prod(data.Nd));
        end

        if VERBOSE, fprintf('FSEG: using %d demodulations (df = %f). Interpolation: %s\n',...
                ncount,data.df,Interpolation);
        end

    % ---------------------------------------------------------
    case 'linear' % Computing all demodulations and interpolating with triang
        % (sometimes is as good, even better, than MFI)

        pmin = min(data.p(:)); pmax = max(data.p(:));
        pixmin = floor(pmin/data.df); pixmax = ceil(pmax/data.df); % min and max index
        pindices = pixmin:pixmax;
        pp = pindices*data.df;      % frequencies in Hz
        Nps = length(pindices); % number of bins

        ncount = 0;
        ims = squeeze(zeros([prod(Nd) Nps]));
        weights = squeeze(zeros([prod(Nd) Nps]));
        if data.Uniform
            for ii = 1:Nps
                temp_im = fftshift(ifftn_n(fftshift((m_hat.*exp(1i*2*pi*t*pp(ii))))));
                ims(:,ii) = temp_im(:);
                weights(:,ii) = triang((data.p(:)-pp(ii))/data.df); % frequency distance to this plane
                ncount = ncount+1;
            end
        else
            for ii = 1:Nps
                temp_im = nufft_adj(m_hat.*exp(1i*2*pi*t*pp(ii)).*data.dcf,data.FTst).';
                ims(:,ii) = temp_im(:);
                weights(:,ii) = triang((data.p(:)-pp(ii))/data.df); % frequency distance to this plane
                ncount = ncount+1;
            end
            ims = ims/sqrt(prod(data.Nd));
        end

        m = reshape(sum(ims.*weights,2),Nd);

        if VERBOSE, fprintf('FSEG: using %d demodulations (df = %f). Interpolation: %s\n',...
                ncount,data.df,Interpolation);
        end

    % ---------------------------------------------------------------
    case 'MFI' % MFI: Computing all demodulations and interpolating with optimal weights

        pmin = min(data.p(:)); pmax = max(data.p(:));
        pixmin = floor(pmin/data.df); pixmax = ceil(pmax/data.df); % min and max index
        pindices = pixmin:pixmax;
        pp = pindices*data.df;      % frequencies in Hz
        Nps = length(pindices); % number of bins

        ims = squeeze(zeros([prod(Nd) Nps]));
        weights = squeeze(zeros([prod(Nd) Nps]));

        % Compute MFI weights in a per pixel base
        pfine = (pixmin:0.25:pixmax)*data.df; % Fine grid of p's
        ppweights = zeros([length(pfine) Nps]);
        A = exp(1i*2*pi*t(:)*pp);
        if ~exist([filename '.mat'],'file') || Nps<=25
            if VERBOSE, fprintf('FSEG: Computing weights.\n');end
            for n = 1:length(pfine)
                b = exp(1i*2*pi*t(:)*pfine(n));
                [bor,~] = lsqr(A,b,[],20);
                ppweights(n,:) = bor.';
            end
            if Nps>25 % Saves weights only if they are slow to compute
                if VERBOSE, fprintf('FSEG: Saving pre-computed weights in %s.mat.\n',filename);end
                saved_Nps = Nps; saved_ndf = data.df;
                save(filename,'ppweights','saved_Nps','saved_ndf');
            end
        else
            if VERBOSE, fprintf('FSEG: Loading pre-computed weights from %s.mat.\n',filename);end
            load(filename);
            if Nps ~= saved_Nps || data.df ~= saved_ndf
                error('Must recompute the weights (inconsistency).'); end
        end

        for n=1:length(data.p(:))
            weights(n,:) = ppweights(round((data.p(n)/data.df-pixmin)*4+1),:);
        end

        ncount = 0;
        if data.Uniform
            for ii = 1:Nps
                temp_im = fftshift(ifftn_n(fftshift((m_hat.*exp(1i*2*pi*t*pp(ii))))));
                ims(:,ii) = temp_im(:);
                ncount = ncount+1;
            end
        else
            for ii = 1:Nps
                temp_im = nufft_adj(m_hat.*exp(1i*2*pi*t*pp(ii)).*data.dcf,data.FTst).';
                ims(:,ii) = temp_im(:);
                ncount = ncount+1;
            end
            ims = ims/sqrt(prod(data.Nd));
        end
        m = reshape(sum(ims.*weights,2),Nd);

        if VERBOSE
            fprintf('FSEG: using %d demodulations (df = %f). Interpolation: %s (fine grid size:%d)\n',...
                ncount,data.df,Interpolation,length(pfine));
        end

    otherwise
        error('Interpolation methos: %s not supported.',Interpolation);

end

end
