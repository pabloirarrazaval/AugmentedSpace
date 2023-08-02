% Generates the figures for the Augmented Space off-resonance paper
%
% Pablo Irarrazaval (June 2023)

%% Preamble
global VERBOSE
VERBOSE = true; % Print messages
addpath('Simulation code'); addpath(genpath('Recon code'));

fig_dir = '..\pdfFigures\'; fig_ext = '.pdf';

% What to do
DO_Sim1d = 0; % 1D simulations
DO_Sim2d = 1; % 2D simulations
DO_Phantom = 0; % Phantom
            LoadPhantom = 1; % load from file because it is already computed
DO_Invivo = 0; % In-vivo
            LoadInvivo = 1; % load from file because it is already computed

SAVE_figure = false;

%% 1D simulation
if DO_Sim1d

    fig_base = 'sim1D_';
    fcnt = 0; f = []; figname = [];

    % Defines reconstruction experiment in data structure
    data = [];
    data.Nd = 128;                       % Image dimensions (Nx)
    data.Tmax = 20e-3;                   % ACQ time in seconds
    data.seq = '2DFT';                   % Defines trajectory
    data.pfun_name = 'Half plane';       % Defines field map
    % data.pfun_name = 'Zero';       % Defines field map

    % Defines simulation parameters in simdata structure
    simdata = [];
    simdata.Nd = data.Nd;                % The dimensions again
    simdata.type = 'Sum';                % Type of simulation
    simdata.objfun_name = 'Centered square'; % Object to simulate

    % Generates simulated data, creates m_hat and updates data and simdata
    Simulate;
    % pw = tukeywin(data.Nd,0.5); % Apodization filter (becaue of the low resolution)
    % m_hat = m_hat.*pw;

    % Shows experiment data
    %-----------------------
    % common parameters for plots
    linewidth = 1.5;

    fcnt = fcnt+1; f{fcnt} = figure(100); figname{fcnt} = 'object'; 
    plot(simdata.x,real(simdata.m),simdata.x,imag(simdata.m),'LineWidth',linewidth);
    % title('Object 1D');
    legend('real','imag');legend('boxoff');
    opt = []; opt.xlabel = '$x$'; opt.ylabel = '$m(x)$';
    opt.yticks = {'1'}; opt.xticks = {'-16' '16'};
    pretty_plot(opt)

    fcnt = fcnt+1; f{fcnt} = figure(101); figname{fcnt} = 'fieldmap';
    plot(simdata.x,data.p,'LineWidth',linewidth);
    axis([-data.Nd/2 data.Nd/2 -20 20])
    % title('Field map [Hz]');
    opt = []; opt.xlabel = '$x$'; opt.ylabel = '$p(x)$';
    opt.yticks = {'-15'};
    pretty_plot(opt);

    fcnt = fcnt+1; f{fcnt} = figure(102); figname{fcnt} = 'timemap';
    plot(data.kx,data.t*1000,'LineWidth',linewidth);
    axis([-0.5 0.5 0 22])
    % title('Time map [ms]')
    opt = []; opt.xlabel = '$k_x$'; opt.ylabel = '$t(k_x)$';
    opt.yticks = {'20'};
    pretty_plot(opt);

    fcnt = fcnt+1; f{fcnt} = figure(103); figname{fcnt} = 'signal';
    plot(data.kx,real(m_hat),data.kx,imag(m_hat),'LineWidth',linewidth);
    axis([-0.5 0.5 -1 3])
    % title('Simulated signal');
    legend('real','imag');legend('boxoff');
    opt = []; opt.xlabel = '$k_x$'; opt.ylabel = '$\hat{m}(k_x)$';
    opt.yticks = {'2.5'}; opt.xticks = {'-0.5' '0.5'};
    pretty_plot(opt);


    % Reconstruct with different algorithms
    %---------------------------------------
    % DFT ****
    m_dft = DFTrecon(m_hat,data);

    fcnt = fcnt+1; f{fcnt} = figure(110); figname{fcnt} = 'dft';
    plot(simdata.x,real(m_dft),simdata.x,imag(m_dft),'LineWidth',linewidth)
    % title('DFT recon')
    legend('real','imag');legend('boxoff');
    opt = []; opt.xlabel = '$x$'; opt.ylabel = '$m(x)$';
    opt.yticks = {'1'}; opt.xticks = {'-16' '16'};
    pretty_plot(opt)

    % FSEG ****
    data.df = 1;    % Frequency step in Hz
    m_fseg = FSEGrecon(m_hat,data,'none'); % No need to interpolate because the
    % field map is already segmented

    fcnt = fcnt+1; f{fcnt} = figure(120); figname{fcnt} = 'fseg';
    plot(simdata.x,real(m_fseg),simdata.x,imag(m_fseg),'LineWidth',linewidth)
    % title('FSEG recon')
    legend('real','imag');legend('boxoff');
    opt = []; opt.xlabel = '$x$'; opt.ylabel = '$m(x)$';
    opt.yticks = {'1'}; opt.xticks = {'-16' '16'};
    pretty_plot(opt)

    % CAS ****
    data.df = 1;        % Sampling period for frequency [Hz]
    data.dt = 0.1e-3;   % Sampling period for time [s]
    data.w = 1;         % For Uniform sampling no need to window data
    for iter = 1:5 % Test different number of iterations
        data.maxiter = iter;
        m_casn = CAS_ls(m_hat,data,[],data.maxiter);

        fcnt = fcnt+1; f{fcnt} = figure(130+iter); figname{fcnt} = sprintf('cas%d',iter);
        plot(simdata.x,real(m_casn),simdata.x,imag(m_casn),'LineWidth',linewidth)
        title(sprintf('CAS recon (%d iters)',data.maxiter))
        legend('real','imag');legend('boxoff');
        opt = []; opt.xlabel = '$x$'; opt.ylabel = '$m(x)$';
        opt.yticks = {'1'}; opt.xticks = {'-16' '16'};
        pretty_plot(opt)
    end

    % DAS ****
    data.Nf = 512;      % Dimension of discrete augmented space
    data.df = 3;        % Sampling period for frequency [Hz]
    data.dt = 1/data.df/data.Nf;   % Sampling period for time [s]
    data.w = 1;         % For Uniform sampling no need to window data
    data.maxiter = 5;
    m_das = DAS_ls(m_hat,data,[],data.maxiter);

    fcnt = fcnt+1; f{fcnt} = figure(140);  figname{fcnt} = 'das';
    plot(simdata.x,real(m_das),simdata.x,imag(m_das),'LineWidth',linewidth)
    title(sprintf('DAS recon (%d iters)',data.maxiter))
    legend('real','imag');legend('boxoff');
    opt = []; opt.xlabel = '$x$'; opt.ylabel = '$m(x)$';
    opt.yticks = {'1'}; opt.xticks = {'-16' '16'};
    pretty_plot(opt)

    if SAVE_figure
        for nn = 1:length(f)
            exportgraphics(f{nn},[fig_dir fig_base figname{nn} fig_ext]);
        end
    end

end

%% 2D simulation
if DO_Sim2d

    fig_base = 'sim2D_';
    fcnt = 0; f = []; figname = []; % store fig handles and names
    times = []; timnames = [];      % store execution times and names

    option = 2; % 1: Resolution, Bipolar gauss, EPI 8, T=20ms
                % 2: Tubes, Bipolar gauss, EPI 8, T=30ms
                % 3: Tubes, Gauss, Spiral 4, T=30ms
                % 4: Tubes, Step field map, Spiral 4, T=30ms

    % Defines reconstruction experiment in data structure
    data = [];
    data.Nd = [128 128];                       % Image dimensions (Nx)
    data.shift = data.Nd/2;              % Position of origin
    switch option
        case 1
            data.Tmax = 20e-3;                   % ACQ time in seconds
            data.seq = 'EPI 8';                   % Defines trajectory
            Nshots = 8;
            data.pfun_name = 'Bipolar gauss';       % Defines field map
            simdata.objfun_name = 'Resolution'; % Object to simulate
        case 2
            data.Tmax = 30e-3;                   % ACQ time in seconds
            data.seq = 'EPI 8';                   % Defines trajectory
            Nshots = 8;
            data.pfun_name = 'Bipolar gauss';       % Defines field map
            simdata.objfun_name = 'Tubes'; % Object to simulate
            PrecomputedMFIweights = 'Recon code\MFIPrecomputedWeights_sim2D';
        case 3
            data.Tmax = 30e-3;                   % ACQ time in seconds
            data.seq = 'Spiral 4';                   % Defines trajectory
            Nshots = 4;
            data.pfun_name = 'Gauss';       % Defines field map
            simdata.objfun_name = 'Tubes'; % Object to simulate
        case 4
            data.Tmax = 30e-3;                   % ACQ time in seconds
            data.seq = 'Spiral 4';                   % Defines trajectory
            Nshots = 4;
            data.pfun_name = 'Half plane';       % Defines field map
            simdata.objfun_name = 'Tubes'; % Object to simulate
            PrecomputedMFIweights = 'Recon code\MFIPrecomputedWeights_sim2Da';
    end

    % Defines simulation parameters in simdata structure
    simdata.Nd = data.Nd;                % The dimensions again
    simdata.type = 'Sum';                % Type of simulation
    

    % Generates simulated data, creates m_hat and updates data and simdata
    Simulate;
    Nsamples = data.L/Nshots; % number of samples per shot

    if data.Uniform
        pw = 1; mask = 1; % no pre-windows, no mask
    else
        % Low pass filters raw data before recon
        pw = tukeywin(2*Nsamples,0.75); pw = repmat(pw(Nsamples+1:end),[1 Nshots]);

        % circle mask to limit FOV to a circle
        [xx,yy] = meshgrid(-data.Nd(1)/2:data.Nd(1)/2-1,-data.Nd(2)/2:data.Nd(2)/2-1);
        mask = (xx.^2+yy.^2)<=(data.Nd(1)/2)^2;
    end

    % Shows experiment data
    %-----------------------
    % common parameters for plots
    lims_obj = [0 1.1]; % same scale for all plots of objects
    lims_err = [0 0.3]; % same scale for all plots of errors

    fcnt = fcnt+1; f{fcnt} = figure(200); figname{fcnt} = 'object'; 
    imshow(abs(simdata.m),lims_obj,'Init',300); set(gca,'YDir','Normal');

    fcnt = fcnt+1; f{fcnt} = figure(201); figname{fcnt} = 'fieldmap';
    imshow(data.p,[],'Init',300); set(gca,'YDir','Normal');

    % need to plot differently non-uniform time and raw data
    if data.Uniform
        fcnt = fcnt+1; f{fcnt} = figure(202); figname{fcnt} = 'timemap';
        imshow(data.t*1000,[],'Init',300);set(gca,'YDir','Normal');

        % Warning: check if signal appears transposed
        fcnt = fcnt+1; f{fcnt} = figure(203); figname{fcnt} = 'signal';
        imshow(abs(m_hat),[],'Init',300);set(gca,'YDir','Normal');
    else
        fcnt = fcnt+1; f{fcnt} = figure(202); figname{fcnt} = 'timemap';
        plot3(data.kx,data.ky,data.t);

        fcnt = fcnt+1; f{fcnt} = figure(203); figname{fcnt} = 'signal';
        plot3(data.kx,data.ky,abs(m_hat),'.');
    end

    % Reconstruct with different algorithms
    %---------------------------------------
    % DFT *****************************************************************
    ii = length(times); tic;
    m_dft = mask.*DFTrecon(m_hat.*pw(:),data);
    times(ii+1) = toc; timnames{ii+1} = 'DFT';

    RMSE_dft = norm(m_dft-simdata.m);

    fcnt = fcnt+1; f{fcnt} = figure(210); figname{fcnt} = 'dft';
    imshow(abs(m_dft),lims_obj,'Init',300);set(gca,'YDir','Normal');

    fcnt = fcnt+1; f{fcnt} = figure(211); figname{fcnt} = 'dfterror';
    imshow(abs(m_dft-simdata.m),lims_err,'Init',300);set(gca,'YDir','Normal');

    % FSEG ****************************************************************
    data.df = 3;    % Frequency step in Hz

    ii = length(times); tic;
    m_fseg = mask.*FSEGrecon(m_hat.*pw(:),data,'none');
    times(ii+1) = toc; timnames{ii+1} = 'FSEG';

    RMSE_fseg = norm(m_fseg-simdata.m);

    ii = length(times); tic;
    m_fseg_linear = mask.*FSEGrecon(m_hat.*pw(:),data,'linear');
    times(ii+1) = toc; timnames{ii+1} = 'FSEG_linear';

    RMSE_fseg_linear = norm(m_fseg_linear-simdata.m);

    ii = length(times); tic;
    m_fseg_mfi = mask.*FSEGrecon(m_hat.*pw(:),data,'MFI',PrecomputedMFIweights);
    times(ii+1) = toc; timnames{ii+1} = 'FSEG_mfi';

    RMSE_fseg_mfi = norm(m_fseg_mfi-simdata.m);

    fcnt = fcnt+1; f{fcnt} = figure(220); figname{fcnt} = 'fseg';
    imshow(abs(m_fseg),lims_obj,'Init',300);set(gca,'YDir','Normal');
    fcnt = fcnt+1; f{fcnt} = figure(221); figname{fcnt} = 'fsegerror';
    imshow(abs(m_fseg-simdata.m),lims_err,'Init',300);set(gca,'YDir','Normal');

    fcnt = fcnt+1; f{fcnt} = figure(222); figname{fcnt} = 'fseglinear';
    imshow(abs(m_fseg_linear),lims_obj,'Init',300);set(gca,'YDir','Normal');
    fcnt = fcnt+1; f{fcnt} = figure(223); figname{fcnt} = 'fseglinearerror';
    imshow(abs(m_fseg_linear-simdata.m),lims_err,'Init',300);set(gca,'YDir','Normal');

    fcnt = fcnt+1; f{fcnt} = figure(224); figname{fcnt} = 'fsegmfi';
    imshow(abs(m_fseg_mfi),lims_obj,'Init',300);set(gca,'YDir','Normal');
    fcnt = fcnt+1; f{fcnt} = figure(225); figname{fcnt} = 'fsegmfierror';
    imshow(abs(m_fseg_mfi-simdata.m),lims_err,'Init',300);set(gca,'YDir','Normal');

    % CAS *****************************************************************
    data.df = 3;        % Sampling period for frequency [Hz]
    data.dt = 0.1e-3;   % Sampling period for time [s]
    data.w = 1;         % For Uniform sampling no need to window data
    data.maxiter = 7;

    ii = length(times); tic;
    m_cas = mask.*CAS_ls(m_hat.*pw(:),data,[],data.maxiter);
    times(ii+1) = toc; timnames{ii+1} = 'CAS';

    RMSE_cas = norm(m_cas-simdata.m);

    fcnt = fcnt+1; f{fcnt} = figure(230); figname{fcnt} = 'cas';
    imshow(abs(m_cas),lims_obj,'Init',300);set(gca,'YDir','Normal');

    fcnt = fcnt+1; f{fcnt} = figure(231); figname{fcnt} = 'caserror';
    imshow(abs(m_cas-simdata.m),lims_err,'Init',300);set(gca,'YDir','Normal');

    % DAS *****************************************************************
    data.Nf = 512;      % Dimension of discrete augmented space
    data.df = 3;        % Sampling period for frequency [Hz]
    data.dt = 1/data.df/data.Nf;   % Sampling period for time [s]
    data.w = 1;         % For Uniform sampling no need to window data
    data.maxiter = 8; 

    ii = length(times); tic;
    m_das = mask.*DAS_ls(m_hat.*pw(:),data,[],data.maxiter);
    times(ii+1) = toc; timnames{ii+1} = 'DAS';

    RMSE_das = norm(m_das-simdata.m);

    fcnt = fcnt+1; f{fcnt} = figure(240); figname{fcnt} = 'das'; 
    imshow(abs(m_das),lims_obj,'Init',300);set(gca,'YDir','Normal');

    fcnt = fcnt+1; f{fcnt} = figure(241); figname{fcnt} = 'daserror';
    imshow(abs(m_das-simdata.m),lims_err,'Init',300);set(gca,'YDir','Normal');

    % Print errors
    fprintf('RMSE DFT = %f FSEG = %f DAS = %f CAS = %f\n',...
        RMSE_dft,RMSE_fseg,RMSE_das,RMSE_cas);

    % Print times
    fprintf('Execution times\n')
    for nn = 1:length(times)
        fprintf('%12s: %5.2f ms\n',timnames{nn},times(nn)*1000);
    end
    fprintf('Execution times\n')
    for nn = 1:length(times)
        fprintf('%12s: %5.2f s\n',timnames{nn},times(nn));
    end


    if SAVE_figure
        for nn = 1:length(f)
            exportgraphics(f{nn},[fig_dir fig_base figname{nn} fig_ext]);
        end
    end

end

%% Phantom acquisition
if DO_Phantom

    fig_base = 'phantom_';
    fcnt = 0; f = []; figname = []; % store fig handles and names
    times = []; timnames = [];      % store execution times and names

    % Speeds things up reconstructing only one coil
    COIL = 0; % 0 for all coils
    FIGURE_TYPE = 'abs'; % 'abs' or 'complex'


    % Defines reconstruction experiment in data structure
    data = [];

    % Reads the data. It contains:
    % - M                           % the raw data
    % - Nd                          % object size
    % - Ncoils, Nsamples, Nshots    % other size variables
    % - Shift                       % to center reconstruction
    % - csm                         % coil sensitivity maps
    % - dcf                         % density compensation function
    % - fm                          % field map
    % - tm                          % time map
    % - k_spx, k_spy                % k-space trajectory
    load('Data\phantom_data.mat');

    % Fills in strcture from loaded variables
    data.Nd = Nd;           % Dimensions (reverse order: z,y,x)
    data.shift = Shift;     % Moves origin
    data.Nc = Ncoils;       % Number of coils
    data.Nsamples = Nsamples; % Number of samples per shot
    data.Nshots = Nshots;   % Number of shots
    data.Uniform = false;               % Uniform or nonuniform sampling
    data.kx = k_spx(:);                 % kx
    data.ky = k_spy(:);                 % ky
    data.dcf = dcf(:);                  % k-space density compensation (ony non-uniform)
    data.L = length(data.kx);           % Number of samples for non-uniform
    data.t = tm(:);                     % Time map
    data.p = fm;                        % Field map
    data.FTst = prepares_nufft(data);   % Prepares NUFFT

    if LoadPhantom % already computed and saved to this file
        load('phantom_figs.mat');
    else
        USEpar = true; % Enables parellel for's in recons
        if USEpar && (COIL==0)
            poolobj = gcp; % Get the current parallel pool and creates it if not available
            NumWorkers = poolobj.NumWorkers;
        else, NumWorkers = 0; end

        % Low pass filters raw data before recon
        pw = tukeywin(2*Nsamples,0.75); pw = repmat(pw(Nsamples+1:end),[1 Nshots]);

        % circle mask to limit FOC to a circle
        [xx,yy] = meshgrid(-Nd(1)/2:Nd(1)/2-1,-Nd(2)/2:Nd(2)/2-1);
        mask = (xx.^2+yy.^2)<=(Nd(1)/2)^2;

        % Function to combine coils using CSM
        if COIL==0 % All coils
            csm_sq = max(0.5,sum(abs(csm).^2,numel(data.Nd)+1)); % 0.5 avoids div by zero
            csmof = @(x) (sum(conj(csm).*x,numel(data.Nd)+1)./csm_sq);
        else % One coil
            csm_sq = (max(0.5,abs(csm(:,:,COIL)).^2)); % 0.5 avoids div by zero
            csmof = @(x) ((conj(csm(:,:,COIL)).*x))./csm_sq;
        end
    end

    % Shows experiment data
    %-----------------------
    % common parameters for plots

    fcnt = fcnt+1; f{fcnt} = figure(301); figname{fcnt} = 'fieldmap'; 
    imshow(data.p,[],'Init',300); set(gca,'YDir','Normal');
%     title('Field map [Hz]');

    fcnt = fcnt+1; f{fcnt} = figure(302); figname{fcnt} = 'timemap'; 
    plot3(data.kx,data.ky,data.t);
%     title('Time map [ms]');


    % Reconstruct with different algorithms
    %---------------------------------------

    % DFT *****************************************************************
    if ~LoadPhantom
        if COIL>0
            s = M(:,:,COIL).*pw;
            m_dft = mask.*nufft_adj(s(:).*data.dcf,data.FTst).'/sqrt(prod(data.Nd));
        else
            m_dft = zeros([data.Nd data.Nc]);
            ii = length(times); tic;
            parfor (n = 1:data.Nc, NumWorkers) % Recon coil by coil
                s = M(:,:,n).*pw;
                m_dft(:,:,n) = mask.*nufft_adj(s(:).*data.dcf,data.FTst).'/sqrt(prod(data.Nd));
            end
            times(ii+1) = toc; timnames{ii+1} = 'DFT';
        end
        % Combine coils
        m_dft_csm = csmof(m_dft);
    end

    fcnt = fcnt+1; f{fcnt} = figure(310); figname{fcnt} = 'dft'; 
    titname = 'DFT recon';
    switch FIGURE_TYPE
        case 'abs'
            imshow(abs(m_dft_csm),[],'Init',300);
            set(gca,'YDir','Normal'); % title(titname);
        case 'complex'
            mkimshowsubplot('complex',m_dft_csm,'common',titname);
    end

    % FSEG NONE ***********************************************************
    if ~LoadPhantom
        data.df = 10;    % Frequency step in Hz
        FSEG_INTERP = 'none';

        if COIL>0
            s =  M(:,:,COIL).*pw;
            m_seg = mask.*FSEGrecon(s(:),data,FSEG_INTERP);
        else
            fprintf('FSEG\n');
            m_seg = zeros([data.Nd data.Nc]);
            ii = length(times); tic;
            parfor (n = 1:data.Nc, NumWorkers) % Recon coil by coil
                fprintf('%d: ',n);
                s =  M(:,:,n).*pw;
                m_seg(:,:,n) = mask.*FSEGrecon(s(:),data,FSEG_INTERP);
            end
            times(ii+1) = toc; timnames{ii+1} = 'FSEG';
            fprintf('\n');
        end
        % Combine coils
        m_seg_csm = csmof(m_seg);
    end

    fcnt = fcnt+1; f{fcnt} = figure(320); figname{fcnt} = 'fseg'; 
    titname = 'FSEG recon';
    switch FIGURE_TYPE
        case 'abs'
            imshow(abs(m_seg_csm),[],'Init',300);
            set(gca,'YDir','Normal'); % title(titname);
        case 'complex'
            mkimshowsubplot('complex',m_seg_csm,'common',titname);
    end

    % FSEG LINEAR ***********************************************************
    if ~LoadPhantom
        data.df = 10;    % Frequency step in Hz
        FSEG_INTERP = 'linear';

        if COIL>0
            s =  M(:,:,COIL).*pw;
            m_seg_linear = mask.*FSEGrecon(s(:),data,FSEG_INTERP);
        else
            fprintf('FSEG linear\n');
            m_seg_linear = zeros([data.Nd data.Nc]);
            ii = length(times); tic;
            parfor (n = 1:data.Nc, NumWorkers) % Recon coil by coil
                fprintf('%d: ',n);
                s =  M(:,:,n).*pw;
                m_seg_linear(:,:,n) = mask.*FSEGrecon(s(:),data,FSEG_INTERP);
            end
            times(ii+1) = toc; timnames{ii+1} = 'FSEG_linear';
            fprintf('\n');
        end
        % Combine coils
        m_seg_linear_csm = csmof(m_seg_linear);
    end

    fcnt = fcnt+1; f{fcnt} = figure(321); figname{fcnt} = 'fseglinear'; 
    titname = 'FSEG linear recon';
    switch FIGURE_TYPE
        case 'abs'
            imshow(abs(m_seg_linear_csm),[],'Init',300);
            set(gca,'YDir','Normal'); % title(titname);
        case 'complex'
            mkimshowsubplot('complex',m_seg_linear_csm,'common',titname);
    end

    % FSEG MFI ***********************************************************
    if ~LoadPhantom
        data.df = 10;    % Frequency step in Hz
        FSEG_INTERP = 'MFI';

        if COIL>0
            s =  M(:,:,COIL).*pw;
            m_seg_mfi = mask.*FSEGrecon(s(:),data,FSEG_INTERP,...
                'Recon code\MFIPrecomputedWeights_phantom');
        else
            fprintf('FSEG mfi\n');
            m_seg_mfi = zeros([data.Nd data.Nc]);
            ii = length(times); tic;
            parfor (n = 1:data.Nc, NumWorkers) % Recon coil by coil
                fprintf('%d: ',n);
                s =  M(:,:,n).*pw;
                m_seg_mfi(:,:,n) = mask.*FSEGrecon(s(:),data,FSEG_INTERP,...
                    'Recon code\MFIPrecomputedWeights_phantom');
            end
            times(ii+1) = toc; timnames{ii+1} = 'FSEG_mfi';
            fprintf('\n');
        end
        % Combine coils
        m_seg_mfi_csm = csmof(m_seg_mfi);
    end

    fcnt = fcnt+1; f{fcnt} = figure(322); figname{fcnt} = 'fsegmfi'; 
    titname = 'FSEG MFI recon';
    switch FIGURE_TYPE
        case 'abs'
            imshow(abs(m_seg_mfi_csm),[],'Init',300);
            set(gca,'YDir','Normal'); % title(titname);
        case 'complex'
            mkimshowsubplot('complex',m_seg_mfi_csm,'common',titname);
    end

    % CAS *****************************************************************
    if ~LoadPhantom
        data.df = 10;       % Sampling period for frequency [Hz]
        data.dt = 0.3e-3;   % Sampling period for time [s]
        data.maxiter = 10;   % Iterations

        % First iteration uses density as weights
        data.w = sqrt(data.dcf);

        if COIL>0
            s =  M(:,:,COIL).*pw;
            m_cas = mask.*CAS_ls(s(:),data,[],1); % First iteration
            data.w = 1;  % For rest of iterations
            s =  M(:,:,COIL).*pw;
            x0 = m_cas;
            m_cas = mask.*CAS_ls(s(:),data,[],data.maxiter-1,[],[],x0(:));
        else
            fprintf('CAS\n');
            m_cas = zeros([data.Nd data.Nc]);
            % First iteration
            fprintf('First iteration with weight of sqrt(dcf):\n')
            ii = length(times); tic;
            parfor (n = 1:data.Nc, NumWorkers) % Recon coil by coil
                fprintf('%d: ',n);
                s =  M(:,:,n).*pw;
                m_cas(:,:,n) = mask.*CAS_ls(s(:),data,[],1);
            end
            fprintf('\n');
            % Rest of iterations
            data.w = 1;
            fprintf('Rest of iterations with weight of one:\n')
            parfor (n = 1:data.Nc, NumWorkers) % Recon coil by coil
                fprintf('%d: ',n);
                s =  M(:,:,n).*pw;
                x0 = m_cas(:,:,n);
                m_cas(:,:,n) = mask.*CAS_ls(s(:),data,[],data.maxiter-1,[],[],x0(:));
            end
            times(ii+1) = toc; timnames{ii+1} = 'CAS';
            fprintf('\n');
        end
        % Combine coils
        m_cas_csm = csmof(m_cas);
    end

    fcnt = fcnt+1; f{fcnt} = figure(330); figname{fcnt} = 'cas'; 
    titname = 'CAS recon';
    switch FIGURE_TYPE
        case 'abs'
            imshow(abs(m_cas_csm),[],'Init',300);
            set(gca,'YDir','Normal'); %title(titname);
        case 'complex'
            mkimshowsubplot('complex',m_cas_csm,'common',titname);
    end

    % DAS *****************************************************************
    if ~LoadPhantom
        data.Nf = 256;                  % Dimension of discrete augmented space
        data.df = 10;                   % Sampling period for frequency [Hz]
        data.dt = 1/data.df/data.Nf;    % Sampling period for time [s]
        data.maxiter = 8;

        % First iteration uses density as weights
        data.w = sqrt(data.dcf);

        if COIL>0
            s =  M(:,:,COIL).*pw;
            m_das = mask.*DAS_ls(s(:),data,[],1); % First iteration
            data.w = 1;  % For rest of iterations
            s =  M(:,:,COIL).*pw;
            x0 = m_das;
            m_das = mask.*DAS_ls(s(:),data,[],data.maxiter-1,[],[],x0(:));
        else
            fprintf('DAS\n');
            m_das = zeros([data.Nd data.Nc]);
            % First iteration
            fprintf('First iteration with weight of sqrt(dcf):\n')
            ii = length(times); tic;
            parfor (n = 1:data.Nc, NumWorkers) % Recon coil by coil
                fprintf('%d: ',n);
                s =  M(:,:,n).*pw;
                m_das(:,:,n) = mask.*DAS_ls(s(:),data,[],1);
            end
            fprintf('\n');
            % Rest of iterations
            data.w = 1;
            fprintf('Rest of iterations with weight of one:\n')
            parfor (n = 1:data.Nc, NumWorkers) % Recon coil by coil
                fprintf('%d: ',n);
                s =  M(:,:,n).*pw;
                x0 = m_das(:,:,n);
                m_das(:,:,n) = mask.*DAS_ls(s(:),data,[],data.maxiter-1,[],[],x0(:));
            end
            times(ii+1) = toc; timnames{ii+1} = 'DAS';
            fprintf('\n');
        end
        % Combine coils
        m_das_csm = csmof(m_das);
    end

    fcnt = fcnt+1; f{fcnt} = figure(340); figname{fcnt} = 'das'; 
    titname = 'DAS recon';
    switch FIGURE_TYPE
        case 'abs'
            imshow(abs(m_das_csm),[],'Init',300);
            set(gca,'YDir','Normal'); % title(titname);
        case 'complex'
            mkimshowsubplot('complex',m_das_csm,'common',titname);
    end

    %------------------------
    load('Data\phantom_reference.mat');
    fcnt = fcnt+1; f{fcnt} = figure(350); figname{fcnt} = 'ref'; 
    titname = 'Reference';
    switch FIGURE_TYPE
        case 'abs'
            imshow(abs(phantom_reference),[],'Init',300);
            set(gca,'YDir','Normal'); % title(titname);
        case 'complex'
            mkimshowsubplot('complex',phantom_reference,'common',titname);
    end

%     if ~LoadPhantom
%         save('phantom_figs.mat','m_*');
%     end

    % Print times
    fprintf('Execution times\n')
    for nn = 1:length(times)
        fprintf('%12s: %5.2f ms\n',timnames{nn},times(nn)*1000);
    end
    fprintf('Execution times\n')
    for nn = 1:length(times)
        fprintf('%12s: %5.2f s\n',timnames{nn},times(nn));
    end

    if SAVE_figure
        for nn = 1:length(f)
            exportgraphics(f{nn},[fig_dir fig_base figname{nn} fig_ext]);
        end
    end

end

%% In-vivo acquisition
if DO_Invivo

    FIGURE_TYPE = 'abszoom'; % 'abs', 'abszoom' or 'complex'
    fig_base = 'invivozoom_'; % change for different types of figures
    fcnt = 0; f = []; figname = []; % store fig handles and names
    times = []; timnames = [];      % store execution times and names

    % Speeds things up reconstructing only one coil
    COIL = 0; % 0 for all coils

    % Defines reconstruction experiment in data structure
    data = [];

    % Reads the data. It contains:
    % - M                           % the raw data
    % - Nd                          % object size
    % - Ncoils, Nsamples, Nshots    % other size variables
    % - Shift                       % to center reconstruction
    % - csm                         % coil sensitivity maps
    % - dcf                         % density compensation function
    % - fm                          % field map
    % - tm                          % time map
    % - k_spx, k_spy                % k-space trajectory
    load('Data\invivo_data.mat');

    % Fills in strcture from loaded variables
    data.Nd = Nd;           % Dimensions (reverse order: z,y,x)
    data.shift = Shift;     % Moves origin
    data.Nc = Ncoils;       % Number of coils
    data.Nsamples = Nsamples; % Number of samples per shot
    data.Nshots = Nshots;   % Number of shots
    data.Uniform = false;               % Uniform or nonuniform sampling
    data.kx = k_spx(:);                 % kx
    data.ky = k_spy(:);                 % ky
    data.dcf = dcf(:);                  % k-space density compensation (ony non-uniform)
    data.L = length(data.kx);           % Number of samples for non-uniform
    data.t = tm(:);                     % Time map
    data.p = fm;                        % Field map
    data.FTst = prepares_nufft(data);   % Prepares NUFFT

    if LoadInvivo % already computed and saved to this file
        load('invivo_figs.mat');
    else
        USEpar = true; % Enables parellel for's in recons
        if USEpar && (COIL==0)
            poolobj = gcp; % Get the current parallel pool and creates it if not available
            NumWorkers = poolobj.NumWorkers;
        else, NumWorkers = 0; end

        % Low pass filters raw data before recon
        pw = tukeywin(2*Nsamples,0.75); pw = repmat(pw(Nsamples+1:end),[1 Nshots]);

        % circle mask to limit FOC to a circle
        [xx,yy] = meshgrid(-Nd(1)/2:Nd(1)/2-1,-Nd(2)/2:Nd(2)/2-1);
        mask = (xx.^2+yy.^2)<=(Nd(1)/2)^2;

        % Function to combine coils using CSM
        if COIL==0 % All coils
            csm_sq = max(0.5,sum(abs(csm).^2,numel(data.Nd)+1)); % 0.5 avoids div by zero
            csmof = @(x) (sum(conj(csm).*x,numel(data.Nd)+1)./csm_sq);
        else % One coil
            csm_sq = (max(0.5,abs(csm(:,:,COIL)).^2)); % 0.5 avoids div by zero
            csmof = @(x) ((conj(csm(:,:,COIL)).*x))./csm_sq;
        end
    end

    % Shows experiment data
    %-----------------------
    % common parameters for plots

    fcnt = fcnt+1; f{fcnt} = figure(401); figname{fcnt} = 'fieldmap'; 
    imshow(data.p,[],'Init',300); set(gca,'YDir','Normal');
%     title('Field map [Hz]');

    fcnt = fcnt+1; f{fcnt} = figure(402); figname{fcnt} = 'timemap'; 
    plot3(data.kx,data.ky,data.t);
%     title('Time map [ms]');


    % Reconstruct with different algorithms
    %---------------------------------------

    % DFT *****************************************************************
    if ~LoadInvivo
        if COIL>0
            s = M(:,:,COIL).*pw;
            m_dft = mask.*nufft_adj(s(:).*data.dcf,data.FTst).'/sqrt(prod(data.Nd));
        else
            m_dft = zeros([data.Nd data.Nc]);
            ii = length(times); tic;
            parfor (n = 1:data.Nc, NumWorkers) % Recon coil by coil
                s = M(:,:,n).*pw;
                m_dft(:,:,n) = mask.*nufft_adj(s(:).*data.dcf,data.FTst).'/sqrt(prod(data.Nd));
            end
            times(ii+1) = toc; timnames{ii+1} = 'DFT';
        end
        % Combine coils
        m_dft_csm = csmof(m_dft);
    end

    fcnt = fcnt+1; f{fcnt} = figure(410); figname{fcnt} = 'dft'; 
    titname = 'DFT recon';
    switch FIGURE_TYPE
        case 'abs'
            imshow(abs(m_dft_csm),[],'Init',300);
            set(gca,'YDir','Normal'); %title(titname);
        case 'abszoom'
            imshowzoomed(abs(m_dft_csm),[42 59 30 30]);
        case 'complex'
            mkimshowsubplot('complex',m_dft_csm,'common',titname);
    end

    % FSEG NONE ****************************************************************
    if ~LoadInvivo
        data.df = 10;    % Frequency step in Hz
        FSEG_INTERP = 'none';

        if COIL>0
            s =  M(:,:,COIL).*pw;
            m_seg = mask.*FSEGrecon(s(:),data,FSEG_INTERP);
        else
            fprintf('FSEG\n');
            m_seg = zeros([data.Nd data.Nc]);
            ii = length(times); tic;
            parfor (n = 1:data.Nc, NumWorkers) % Recon coil by coil
                fprintf('%d: ',n);
                s =  M(:,:,n).*pw;
                m_seg(:,:,n) = mask.*FSEGrecon(s(:),data,FSEG_INTERP);
            end
            times(ii+1) = toc; timnames{ii+1} = 'FSEG';
            fprintf('\n');
        end
        % Combine coils
        m_seg_csm = csmof(m_seg);
    end

    fcnt = fcnt+1; f{fcnt} = figure(420); figname{fcnt} = 'fseg'; 
    titname = 'FSEG recon';
    switch FIGURE_TYPE
        case 'abs'
            imshow(abs(m_seg_csm),[],'Init',300);
            set(gca,'YDir','Normal'); %title(titname);
        case 'abszoom'
            imshowzoomed(abs(m_seg_csm),[42 59 30 30]);
        case 'complex'
            mkimshowsubplot('complex',m_seg_csm,'common',titname);
    end

    % FSEG LINEAR ****************************************************************
    if ~LoadInvivo
        data.df = 10;    % Frequency step in Hz
        FSEG_INTERP = 'linear';

        if COIL>0
            s =  M(:,:,COIL).*pw;
            m_seg_linear = mask.*FSEGrecon(s(:),data,FSEG_INTERP);
        else
            fprintf('FSEG\n');
            m_seg_linear = zeros([data.Nd data.Nc]);
            ii = length(times); tic;
            parfor (n = 1:data.Nc, NumWorkers) % Recon coil by coil
                fprintf('%d: ',n);
                s =  M(:,:,n).*pw;
                m_seg_linear(:,:,n) = mask.*FSEGrecon(s(:),data,FSEG_INTERP);
            end
            times(ii+1) = toc; timnames{ii+1} = 'FSEG_linear';
            fprintf('\n');
        end
        % Combine coils
        m_seg_linear_csm = csmof(m_seg_linear);
    end

    fcnt = fcnt+1; f{fcnt} = figure(421); figname{fcnt} = 'fseglinear'; 
    titname = 'FSEG linear recon';
    switch FIGURE_TYPE
        case 'abs'
            imshow(abs(m_seg_linear_csm),[],'Init',300);
            set(gca,'YDir','Normal'); %title(titname);
        case 'abszoom'
            imshowzoomed(abs(m_seg_linear_csm),[42 59 30 30]);
        case 'complex'
            mkimshowsubplot('complex',m_seg_linear_csm,'common',titname);
    end

    % FSEG MFI ****************************************************************
    if ~LoadInvivo
        data.df = 10;    % Frequency step in Hz
        FSEG_INTERP = 'MFI';

        if COIL>0
            s =  M(:,:,COIL).*pw;
            m_seg_mfi = mask.*FSEGrecon(s(:),data,FSEG_INTERP,...
                'Recon code\MFIPrecomputedWeights_invivo');
        else
            fprintf('FSEG\n');
            m_seg_mfi = zeros([data.Nd data.Nc]);
            ii = length(times); tic;
            parfor (n = 1:data.Nc, NumWorkers) % Recon coil by coil
                fprintf('%d: ',n);
                s =  M(:,:,n).*pw;
                m_seg_mfi(:,:,n) = mask.*FSEGrecon(s(:),data,FSEG_INTERP,...
                    'Recon code\MFIPrecomputedWeights_invivo');
            end
            times(ii+1) = toc; timnames{ii+1} = 'FSEG_mfi';
            fprintf('\n');
        end
        % Combine coils
        m_seg_mfi_csm = csmof(m_seg_mfi);
    end

    fcnt = fcnt+1; f{fcnt} = figure(422); figname{fcnt} = 'fsegmfi'; 
    titname = 'FSEG MFI recon';
    switch FIGURE_TYPE
        case 'abs'
            imshow(abs(m_seg_mfi_csm),[],'Init',300);
            set(gca,'YDir','Normal'); %title(titname);
        case 'abszoom'
            imshowzoomed(abs(m_seg_mfi_csm),[42 59 30 30]);
        case 'complex'
            mkimshowsubplot('complex',m_seg_mfi_csm,'common',titname);
    end

    % CAS *****************************************************************
    if ~LoadInvivo
        data.df = 10;       % Sampling period for frequency [Hz]
        data.dt = 0.3e-3;   % Sampling period for time [s]
        data.maxiter = 8;   % Iterations

        % First iteration uses density as weights
        data.w = sqrt(data.dcf);

        if COIL>0
            s =  M(:,:,COIL).*pw;
            m_cas = mask.*CAS_ls(s(:),data,[],1); % First iteration
            data.w = 1;  % For rest of iterations
            s =  M(:,:,COIL).*pw;
            x0 = m_cas;
            m_cas = mask.*CAS_ls(s(:),data,[],data.maxiter-1,[],[],x0(:));
        else
            fprintf('CAS\n');
            m_cas = zeros([data.Nd data.Nc]);
            % First iteration
            fprintf('First iteration with weight of sqrt(dcf):\n')
            ii = length(times); tic;
            parfor (n = 1:data.Nc, NumWorkers) % Recon coil by coil
                fprintf('%d: ',n);
                s =  M(:,:,n).*pw;
                m_cas(:,:,n) = mask.*CAS_ls(s(:),data,[],1);
            end
            fprintf('\n');
            % Rest of iterations
            data.w = 1;
            fprintf('Rest of iterations with weight of one:\n')
            parfor (n = 1:data.Nc, NumWorkers) % Recon coil by coil
                fprintf('%d: ',n);
                s =  M(:,:,n).*pw;
                x0 = m_cas(:,:,n);
                m_cas(:,:,n) = mask.*CAS_ls(s(:),data,[],data.maxiter-1,[],[],x0(:));
            end
            times(ii+1) = toc; timnames{ii+1} = 'CAS';
            fprintf('\n');
        end
        % Combine coils
        m_cas_csm = csmof(m_cas);
    end

    fcnt = fcnt+1; f{fcnt} = figure(430); figname{fcnt} = 'cas'; 
    titname = 'CAS recon';
    switch FIGURE_TYPE
        case 'abs'
            imshow(abs(m_cas_csm),[],'Init',300);
            set(gca,'YDir','Normal'); %title(titname);
        case 'abszoom'
            imshowzoomed(abs(m_cas_csm),[42 59 30 30]);
        case 'complex'
            mkimshowsubplot('complex',m_cas_csm,'common',titname);
    end

    % DAS *****************************************************************
    if ~LoadInvivo
        data.Nf = 256;                  % Dimension of discrete augmented space
        data.df = 10;                   % Sampling period for frequency [Hz]
        data.dt = 1/data.df/data.Nf;    % Sampling period for time [s]
        data.maxiter = 8;

        % First iteration uses density as weights
        data.w = sqrt(data.dcf);

        if COIL>0
            s =  M(:,:,COIL).*pw;
            m_das = mask.*DAS_ls(s(:),data,[],1); % First iteration
            data.w = 1;  % For rest of iterations
            s =  M(:,:,COIL).*pw;
            x0 = m_das;
            m_das = mask.*DAS_ls(s(:),data,[],data.maxiter-1,[],[],x0(:));
        else
            fprintf('DAS\n');
            m_das = zeros([data.Nd data.Nc]);
            % First iteration
            fprintf('First iteration with weight of sqrt(dcf):\n')
            ii = length(times); tic;
            parfor (n = 1:data.Nc, NumWorkers) % Recon coil by coil
                fprintf('%d: ',n);
                s =  M(:,:,n).*pw;
                m_das(:,:,n) = mask.*DAS_ls(s(:),data,[],1);
            end
            fprintf('\n');
            % Rest of iterations
            data.w = 1;
            fprintf('Rest of iterations with weight of one:\n')
            parfor (n = 1:data.Nc, NumWorkers) % Recon coil by coil
                fprintf('%d: ',n);
                s =  M(:,:,n).*pw;
                x0 = m_das(:,:,n);
                m_das(:,:,n) = mask.*DAS_ls(s(:),data,[],data.maxiter-1,[],[],x0(:));
            end
            times(ii+1) = toc; timnames{ii+1} = 'DAS';
            fprintf('\n');
        end
        % Combine coils
        m_das_csm = csmof(m_das);
    end

    fcnt = fcnt+1; f{fcnt} = figure(440); figname{fcnt} = 'das'; 
    titname = 'DAS recon';
    switch FIGURE_TYPE
        case 'abs'
            imshow(abs(m_das_csm),[],'Init',300);
            set(gca,'YDir','Normal'); %title(titname);
        case 'abszoom'
            imshowzoomed(abs(m_das_csm),[42 59 30 30]);
        case 'complex'
            mkimshowsubplot('complex',m_das_csm,'common',titname);
    end


    %------------------------
    load('Data\invivo_reference.mat');
    fcnt = fcnt+1; f{fcnt} = figure(450); figname{fcnt} = 'ref'; 
    titname = 'Reference';
    switch FIGURE_TYPE
        case 'abs'
            imshow(abs(invivo_reference),[],'Init',300);
            set(gca,'YDir','Normal'); %title(titname);
        case 'abszoom'
            imshowzoomed(abs(invivo_reference),[42 59 30 30]);
        case 'complex'
            mkimshowsubplot('complex',invivo_reference,'common',titname);
    end

%     if ~LoadInvivo
%         save('invivo_figs.mat','m_*');
%     end

    % Print times
    fprintf('Execution times\n')
    for nn = 1:length(times)
        fprintf('%12s: %5.2f ms\n',timnames{nn},times(nn)*1000);
    end
    fprintf('Execution times\n')
    for nn = 1:length(times)
        fprintf('%12s: %5.2f s\n',timnames{nn},times(nn));
    end

    if SAVE_figure
        for nn = 1:length(f)
            exportgraphics(f{nn},[fig_dir fig_base figname{nn} fig_ext]);
        end
    end

end
