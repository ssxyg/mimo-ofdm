clc;clear;close all;
% set parameters
params.M = 4; %1:BPSK, 2:QPSK, 4:16QAM, 6: 64QAM
params.R = '1/2'; % 1/3 1/2 3/4
params.N = 2048;
params.CP = 152;
params.Nsym = 12; % 60 data symbols per frame
params.Nsubc = 1408; %352 data sub carriers per symbol
params.pilotSpace = 5;  % p d d d d p d d d d ...
params.windowLen = 18;  % length for window process
params.hasSignal = true;  % support dynamic mcs decoding or not
params.agcLen = 16*100;  % length of agc seq, must be multiple of 16
params.pilotFactor =1i;
params.tr_enable = false; % lost 1 db
params.tr_L = 4;
params.tr_A = 0.65;
params.tr_subc = 64;
params.tr_clip = false;
params.tr_clip_factor = 0.75; % 2.75;

params.sim_fs = 20e6;
params.sim_fp = true;

params.plot_spectrum = true;
params.plot_iq = true;
params.plot_ccdf = true;

% prepare encoder for reuse
cnfg.name='LDPC Encoder';
cnfg.operation = 1;
cnfg.standard = 1;
cnfg.algorithm = 1;
cnfg.packing = 1;
cnfg.output_parity_check = 0;
cnfg.bypass = 0;
cnfg.ip_quant_mode =0;

addpath('./ldpc_v2_0_bitacc_cmodel_nt64');
ldpcEnc = ldpc_v2_0_bitacc(cnfg);
rng(11);

% TxGenChanBits
[chanBits,srcBytes,txBlkNum,txSize] = TxGenChanBits(params,ldpcEnc);


% OFDM data modulation
[dataSyms,pilots,refTxQam] = TxGenDataSyms(chanBits,params);

if params.hasSignal
    userId = 255; % 255 means beacon from AP
    N = 50;
    [sigSym, sigByte] = TxGenSigSym(N,userId,params,ldpcEnc);
    dataSyms = [sigSym,dataSyms];
end

% OFDM CP and Window
frameDataPart = TxGenFrameData(dataSyms,params);

% Add preamble
framePreamble = TxGenFramePre(params);

txFrame = [framePreamble;frameDataPart];
save('TX_out.mat');

tx_iq_data = fi(txFrame,1,12,11);
fid = fopen('tx_iq.txt.','w');
for ii = 1:length(tx_iq_data)
    fprintf(fid,'%s%s\n',hex(imag(tx_iq_data(ii))),hex(real(tx_iq_data(ii))));
end
fclose(fid);

if params.plot_iq
    figure;
    subplot(311);
    plot(real(txFrame))
    subplot(312);
    plot(imag(txFrame))
    subplot(313);
    plot(abs(txFrame).^2)
end


if params.plot_spectrum
    hsa = dsp.SpectrumAnalyzer('SampleRate',params.sim_fs);
    step(hsa,[txFrame]);
%     hsa(txFrame);
elseif params.plot_ccdf
    hsa = dsp.SpectrumAnalyzer('SampleRate',params.sim_fs);
    hsa.CCDFMeasurements.Enable = true;
    hsa.CCDFMeasurements.PlotGaussianReference = true;
    step(hsa,frameDataPart);
end
n = length(txFrame);  
power_fre=abs(fft(txFrame))./(params.sim_fs/n);
powerdb_fre=10*log10(power_fre);
                       
fshift = (-n/2:n/2-1)*(params.sim_fs/n);
 yshift = fftshift(powerdb_fre);
figure;
plot(fshift,yshift);
