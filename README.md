# Scripts for running DAS and CAS reconstruction
Pablo Irarrazaval, August 2023

## 1 Introduction
This folder contains all Matlab functions and scripts to run the reconstructions. The recons themselve are stored in folder `Recon code`. The four reconstructions are 
- `DFTrecon.m`: no off-resonance correction. It calls FFT directly.
- `FSEGrecon.m`: Frequency Segmented reconstruction, including interpolated versions.
- `DAS_ls.m`: Discrete Augmented Space reconstruction.
- `CAS_ls.m`: Continuous Augmented Space reconstruction.
  
A good starting point is to read and run the script `generate_figures.m`. This script generates all the figures of the paper. After setting up the paths, there are several boolean variables to help select only parts of the script. `DO_Sim1d`, `DO_Sim2d`, `DO_Phantom`, and `DO_Invivo` choose which simulation to make. All of them will reconstruct the data using the four recons, and all variants of Frequency Segmented. It probably is easier and faster to start with the simulations first. Make sure to set `LoadPhantom = 0` and `LoadInvivo = 0` otherwise, it will read the reconstruction from a file instead of doing it.

### Fourier transforms
To compute the Fourier transform in the uniform sampling case, we use the functions `fftn_n()` and `ifft_n()`. They are `fftn()` and `ifftn()` from Matlab but normalized such that both have a factor `1/sqrt(N)` in front. They are both defined in folder `Recon code\Aux scripts`.

For the non-uniform sampling we use the NUFFT library from Jeff Fessler, Brad Sutton, and Yingying Zhang in folder `Recon code\fessler_nufft`. We use `nufft()` for the direct DFT and `nufft_adj` for the adjoint DFT. The inverse DFT is obtained by computing the adjoint of the raw data weighted by the Density Compensation Function.

## 2 The recon scripts
The four recons are called by:

`m = recon(m_hat,data)`

where
- `m` is the reconstructed object.
- `m_hat` is the raw data (acquired). For uniform sampling is a matrix with the size of the problem. For non-uniform, is a matrix (Number of shots times number of samples per shot) or a vector with the samples. Their k-space positions are given in data.
- `data` is a structure with the problem parameters.
  
Fields of the data structure:
- `Nd`: Dimensions of the object (number of elements equal to the number of dimensions). `[Nx]` or `[Ny Nx]` or `[Nz Ny Nx]`. Note that 3D has not been tested.
- `Nf`: The number of samples for the augmented dimension (f and t).
- `Uniform`: true for uniform sampling and false for non-uniform sampling.
- `p`: field map in Hertz, size Nd.
- `t`: time map in seconds. One entry per k-pace sample.
- `df`: sampling period in the augmented frequency dimension
- `dt`: sampling period in the augmented time dimension

Only needed for non-uniform sampling:
- `L`: Number of k-space samples in m_hat.
- `kx` (and `ky`): k-space sample positions, normalized to -0.5 to 0.5.
- `dcf`: Density Compensation Function.
- `w`: Weight function used to compute a weighted least square solution.
- `shift`: Defines the center of object (typically `Nd/2`). 
- `FTst`: Structure required by NUFFT that must be initialized with `data.FTst = prepares_nufft(data)`.

### 2.1 DFT
`m = DFTrecon(m_hat,data)`: The DFT recon simply calls `ifft_n` for uniform sampling and the `nufft_adj` for the non-uniform sampling.

### 2.2 Frequency Segmented
`m = FSEGrecon(m_hat,data,[optional pars])`: Frequency Segmented recon. The default FSEGrecon calculates the Frequency Segmented reconstruction with no interpolation. FSEGrecon can receive an optional parameter to indicate the type of interpolation: `'none'` (default), `'linear'`, `'MFI'` and `'fix'`. MFI is the Multi-Frequency Interpolation. It needs to compute the weights per pixel for the given field map in a fine grid. It uses a 4x grid. Since this computation is slow, it also has the option to read those weights from a file previously stored. The name of the file is a second optional parameter. fix is for debugging. It demodulates to a fix frequency given as a second optional parameter.

### 2.3 Continuous Augmented Space
`[m,varargout] = CAS_ls(m_hat,data,[varargin])`: The CAS_ls recon. It finds the least square optimum to `||Em-s||` using `lsqr` from Matlab. The direct and adjoint operators are defined in `CE` and `CEadj`, which are called by lsqr. CAS_ls can receive optional input or generate optional output parameters. These are passed or received directly from lsqr.

The direct and adjoint operators are called with parameters in the structure reconopts: `y = CE(r,reconopts)` and `y = CEadj(r,reconopts)`. Fiels of `reconopts` are:
- `Nd`, `Uniform`, `rawdata_dim` (and `FTst`, `w` only for Non-uniform)
- For CE: `tix`, `dt`,`p`
- For CEadj: `pix`, `df`, `t`

Note: `mylsqr` is a copy of the original Matlab's `lsqr`. `mylsqr` outputs the cost function value per iteration (it is nice to have some feedback on what is going on while waiting!). Change it back to `lsqr` to have the standard behavior.

### 2.4 Discrete Augmented Space
`[m,varargout] = DAS_ls(m_hat,data,[varargin])`: The DAS_ls recon. It finds the least square optimum to `||Em-s||` using `lsqr` from Matlab. The direct and adjoint operators are defined in `DE` and `DEadj`, which are called by lsqr. DAS_ls can receive optional input or generate optional output parameters. These are passed or received directly from lsqr.

For DAS, one must have `Nf*df*dt = 1`, to fulfill Nyquist criteria.

The direct and adjoint operators are called with parameters in the structure reconopts: `y = DE(r,reconopts)` and `y = DEadj(r,reconopts)`. Fiels of `reconopts` are:
- `Nd`, `Nf`, `Uniform`, `rawdata_dim`, `linpix`, `lintix` (and `FTst`, `w` only for Non-uniform)
- For DE: `idxpix`
- For DEadj: `idxtix`

Note: `mylsqr` is a copy of the original Matlab's `lsqr`. `mylsqr` outputs the cost function value per iteration (it is nice to have some feedback on what is going on while waiting!). Change it back to `lsqr` to have the standard behavior.

## 3 Making indices
Internally, the recon scripts calls the auxiliary function `pixtix`.

`data = pixtix(data,'op')` is used to convert real-valued frequencies from the field map (p) or time from the time map (t) to indices into matrices (`pix` and `tix`). The operation can be any of the following. Natural order refers to indices (N of them) with the origin in the middle.
- `'create +'`:   Creates `pix` and `tix` in natural order from 1 to N.
- `'create -'`:   Creates `pix` and `tix` in natural order (simple discretization, it includes negatives)
- `'shift x'`:    fftshift's in the image/k-space direction already existing `pix` and `tix`.
- `'shift f'`:    fftshift's in the frequency/time direction already existing `pix` and `tix`.
- `'linear'`:     Creates `linpix` and `lintix`. Linearized indices to arrays to speed up access.
- `'avail'`:      Creates `idxpix` and `idxtix`. Indices of `pix` and `tix` with values different from zero. It is used to avoid computing fft's for only zeros.
