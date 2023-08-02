function m = DFTrecon(m_hat,data)
% Simple wrapper for calling the DFT on m_hat.

if data.Uniform
    m = fftshift(ifftn_n(fftshift(m_hat)));
else
    m = nufft_adj(m_hat.*data.dcf,data.FTst).'/sqrt(prod(data.Nd));
end