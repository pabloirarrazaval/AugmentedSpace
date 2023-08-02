% Script that calls all the simulation procedures, to be called from
% generate_figures
%
% Pablo Irarrazaval (June 2023)

% Generates the k-space trajectory and time map
data = mk_traj(data);
% Generates the function with the test object
simdata = mk_objfun(simdata);
% Generates the field map
data = mk_pfun(simdata,data); % Needs simdata for simdata.x

if ~data.Uniform, data.FTst = prepares_nufft(data); end % Prepares NUFFT

% Simulates and stores the acquisition (if not already pre-computed)
if length(simdata.Nd)==1 % If 1D do it right away, otherwise check if it is precomputed
    name = sprintf('%s, %.2fms, Seq %s, %s',...
        simdata.objfun_name,1000*data.Tmax,data.seq,data.pfun_name);
    m_hat = ACQsimulate(data,simdata,name);
else
    dims = sprintf('%dx%d',data.Nd(1),data.Nd(2));
    filename = sprintf('Pre-computed simulations\\%s, size=%s, Tmax=%.2fms, Seq=%s, fm=%s.mat',simdata.objfun_name, ...
        dims,1000*data.Tmax,data.seq,data.pfun_name);
    if isfile(filename)
        load(filename);
        fprintf('Loaded %s\n',filename);
    else
        name = sprintf('%s, %.2fms, Seq %s, %s',simdata.objfun_name,1000*data.Tmax,data.seq,data.pfun_name);
        m_hat = ACQsimulate(data,simdata,name); save(filename,'m_hat');
        fprintf('Saved %s\n',filename);
    end
end

