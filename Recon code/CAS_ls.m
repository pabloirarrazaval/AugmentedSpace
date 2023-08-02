function [varargout] = CAS_ls(m_hat,data,varargin)
% Continuous Augmented Space least squares
%
% Off-Resonance aware reconstruction. It augments the dimensions adding
% time frequency to image space and time (or kf) to Fourier space, 
% continuously. Then calls lsqr with an efficiente implementation of the matrix
% multiplication.
%
% CAS_ls(m_hat,data,varargin) solves m_hat = A*m
%       m_hat:     acquired data (hat because it is in k_space)
%       data:      information about the problem
%                   .Nd: dimensions
%                   .Nf: size of augmentation
%                   .L:  Number of k-space samples
%                   .FTst: structure with information for NUFFT (only used for non-uniform)
%       varargin:  inputs that will be passed directly to lsqr
%
%       varargout: output (same as lsqr)

%
% Pablo Irarrázaval
% 17 jan 2023: Based on DAS_ls, but computing the adjoint using Conjugate
% Phase
% 1 mar 2023: Uses external definition of operators (CE and CEadj)

global VERBOSE

% Starts by getting number of dimensions and size
Nd = data.Nd;
DIM = length(Nd);
if (~(DIM==2) && ~data.Uniform), error('Non-uniform only works in 2D.'); end
if DIM == 1, Nd = [Nd 1]; end  % to use reshape command (it needs at least two components)
if ~data.Uniform, L = data.L; end

% Check that important information is stored in data
if isempty(data.Uniform) || isempty(data.t)
    error('Trajectory not defined (empty data.t).'); end
if isempty(data.p)
    error('Field map not defined (empty data.p).'); end
if ~data.Uniform && isempty(data.FTst)
    error('data.FTst needed for Non-uniform reconstruction.'); end

% Discretize (finely to simulate continuity) the field and time map
data = pixtix(data,'create -');

% Fills in recon structure
reconopts.Nd = Nd;
reconopts.Uniform = data.Uniform;
if data.Uniform
    reconopts.rawdata_dim = Nd;
else
    reconopts.rawdata_dim = [L 1];
    reconopts.FTst = data.FTst;
end
reconopts.tix = data.tix;
reconopts.pix = data.pix;
reconopts.dt = data.dt;
reconopts.df = data.df;
reconopts.p = data.p;
reconopts.t = data.t;
reconopts.w = data.w;

if VERBOSE
    afuncalls = 0; % records number of callings to afun
    if data.Uniform, funname = 'FFT'; else funname = 'NUFFT'; end
    fprintf('CAS_ls: using df = %f Hz and dt = %f ms.\n',data.df,data.dt*1000);
end

if data.Uniform
    b = fftshift(m_hat); % more efficiente to work in FFT order in optimization
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % The optimization min ||Em - m_hat||_2
    [varargout{1:nargout}] = lsqr(@afun,b(:),varargin{:});
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    varargout{1} = fftshift(reshape(varargout{1},Nd)); % Go back to Natural order
else
    b = m_hat(:).*data.w;
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % The optimization min ||Em - m_hat||_2
    [varargout{1:nargout}] = mylsqr(@afun,b,varargin{:});
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    varargout{1} = reshape(varargout{1},Nd);
end

if VERBOSE, fprintf('\n'); end % end count of callings to afun

%-------------------------------------
    function y = afun(r,flag)

%         if VERBOSE 
%             if afuncalls == 0, afuncalls = 1; fprintf('%s calls: %3d ',funname,afuncalls);
%             else, afuncalls = afuncalls + 1; fprintf('\b\b\b\b%3d ',afuncalls); end
%         end

        if strcmp(flag,'notransp') % Compute A*r

            y = CE(r,reconopts);
     
        elseif strcmp(flag,'transp') % Compute A'*r

            y = CEadj(r,reconopts);

        end

    end


end