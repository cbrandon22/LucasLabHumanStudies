function [wvlt_amp wvlt_phase]  = morletwave(frequencies,c,signal,fs,bl_range,varargin)
% wvlt  = myWaveletTrials_fast4(frequencies,c,signal,fs,bc)
% Compute absolute of Morlet wavelet spectrograms.
% Based on myWaveletTrials_fast5, but with corrected normalization
% 
% Inputs:
% frequencies: a vector of the frequencies for the wavelet transform
% c: the wavelet constant - the ratio between the window width and the cycle length.
% c: the wavelet constant - the ratio between the mean (the center
% frequency) and the standard deviation (the window width) of the frequency
% Gaussian window, or the number of oscillation cycles within a +- 1
% standard deviation of the temporal Gaussian window times pi: 
% c = fc / sig_f = 2pi*sig_t / T.   
%     c = 0 yields the original time-courses (a zero width time window);
%     c -->inf yields the Fourier transform
%     default is 12 (equals 8 in Analyzer). 
% signal: the data matrix - each row is a time course of a single trial or channel.
% fs: Sampling frequency of the signal. Defualt = 1024Hz.
% bl_range: Baseline Correction - 0=no BC; 1 = BC. [a b] = specify range. Default=100:200
% 
% Optional inputs:
% 'waitbar' followed by 'on' or 'off' - controlls whether a waitvar
% appears (default is 'on').
% 'showtext' followed by 'on' or 'off' - controlls whether output text is
% sent to the console. Default is 'off'.
% 'norm' followed by 'on' or 'off' - determines whether Hilbert-compatible
% normalization is used. When set to 'on', the output wavelet amplitudes
% will reflect the actual enevelope amplitude of the band-passed signal, as
% with the Hilbert transform (otherwise, the result is multiples by a
% constant). Recommended to be 'on' but default is 'off' for backward compatibility. 
% 'gpu' followed by 'on' or 'off' - specifies whether to use the GPU for
% the wavelet computation. Testing showed a 20-50% decrease in run-time for
% large inputs, however in certain cases (especially small inputs) run-time
% may actually be slower. Compare the run time for a single run and choose the faster option
% before using in a loop. 
%
% Output:
% wvlt_amp: the time-frequency wavelet matrices for all trials or channels.
%     It is of size (# of frequencies) x (# of time points) x (# of trials or channels)
% wvlt_phase: the time-phase wavelet matrices for all trials or channels.
%     It is the same size and format as wvlt_amp.
% 
% Written by Orr Tomer and Leon Deouell 2007
% Revised by Alon Keren 13.1.2008 (changed to conv2)
% Revised by Edden Gerber 2010-2011
%
% Alon Keren 25.5.08: Changed time window indices to fit short epochs and
% took the baseline correction out of the loop.
% Alon Keren 23.6.08: Added division by 2 in the exponent of the Gaussian
% envelope.
% Alon Keren 1.1.09: Changed display to a waitbar
% Alon Keren 2.1.09: Shortened convolution windows to increase speed
% Alon Keren 30.8.09: Corrected the normalization of the wavelets
% Edden Gerber 17.10.2010: Added an optional input parameter list. One
% input is currently supported: 'waitbar', followed by 'on' or 'off'.
% Default is 'on' for backward compatibility. 
% Edden Gerber 17.10.2010: Changed the input parameter 'bc' to 'bl_range', 
% which will specify a base-line range ( [start end] ). To maintain backward 
% compatibility, if the the value is a boolean, the default [100 200] baseline 
% is used. 
% Edden Gerber 16.2.2011: Added an optional input parameter: 'showtext'
% followed by 'on' or 'off' - controls whether textual output is displayed
% during function processing. Default is 'off'. 
% Edden Gerber 16.2.2011: Added parallel computing option. Enter optional
% parameter 'parallel' followed by 'on'. Default is 'off'. If a
% "matlabpool" is already open, then this function will use it and keep it
% open. If not, this function will open a pool and close it at termination.
% Edden Gerber 2.6.2011: Removed the parallel computing option, as the conv
% function is already configured to run in this mode by default. 
% Edden Gerber 2.6.2011: Added the additional output variable wvlt_phase -
% for analysis of intantaneous phase. 
% Edden Gerber 14.6.2011: If there is only one trial to be processed, the
% waitbar shows progression of wavelet frequencies. 
% Edden Gerber 19.9.2012: Added optional parameters 'norm' - when set to
% 'on', the wavelet will be normalized by the integral of the absolute
% wavelet, leading to an output amplitude equivalet to that of the Hilbert
% transform. 
% Edden Gerber 10.2.2013: Added description of optional input variables in
% function help. 
% Edden Gerber 10.4.2013: Added GPU processing support. 

if mod(size(varargin),2) == 1
    error('Wrong number of input arguments.');
end
if (length(bl_range) ~= 2) && (bl_range ~=0)
    error('bl_range should be 0 or a two-element array');
end

Normalize = false;
ShowWB = true;
ShowText = false;
UseGpu = false;
arg  = 1;
while arg < size(varargin,2)
    switch varargin{arg}
        case 'normalize'
            if strcmp(varargin{arg+1}, 'on')
                Normalize = true;
            elseif strcmp(varargin{arg+1}, 'off')
                Normalize = false;
            else
                error('Illegal parameter value: use ''on'' or ''off''.');
            end
            arg = arg + 2;
        case 'waitbar'
            if strcmp(varargin{arg+1}, 'on')
                ShowWB = true;
            elseif strcmp(varargin{arg+1}, 'off')
                ShowWB = false;
            else
                error('Illegal parameter value: use ''on'' or ''off''.');
            end
            arg = arg + 2;
        case 'showtext'
            if strcmp(varargin{arg+1}, 'on')
                ShowText = true;
            elseif strcmp(varargin{arg+1}, 'off')
                ShowText = false;
            else
                error('Illegal parameter value: use ''on'' or ''off''.');
            end
            arg = arg + 2;
        case 'gpu'
            if strcmp(varargin{arg+1}, 'on')
                UseGpu = true;
            elseif strcmp(varargin{arg+1}, 'off')
                UseGpu = false;
            else
                error('Illegal parameter value: use ''on'' or ''off''.');
            end
                arg = arg + 2;
        otherwise
            error(['Unknown optional argument name: ' varargin{arg} '.']);
    end
end

if ~exist('c','var')
    c = 12; 
end
if ~exist('fs','var')
    fs = 1024; 
end
if islogical(bl_range) || (length(bl_range) == 1)
    if bl_range
        bl_range = [100 200];
    else
        bl_range = 0;
    end
end

[trials, Nts] = size(signal);
Nf = length(frequencies);
    
sig_f = frequencies/c;  
sig_t = 1./(2*pi*sig_f);    % shouldn't 2*pi be sqrt-ed? Alon 12.10.08
twinw = 4 * sig_t'; % Convolution time window half width at each frequency
winw = floor(twinw*fs);
win = winw*ones(1,Nts) + ones(Nf,1)*(1:Nts);
% t=-1:(1/fs):1;
[maxf,imax] = min(frequencies);
t= -twinw(imax):(1/fs):twinw(imax);
Ntw = length(t);    % # of time points in each wavelet
wlim = ceil(Ntw/2) + winw * [-1 1];
%morlet  = (pi*sig_t*sig_t)^(-0.25)*exp(2*i*pi*f*t).*exp(-power(t,2)/(sig_t^2));
% normalize = power(pi*power(sig_t,2),-0.25);   % Alon 23.6.08
% normalize = pi^-0.25 * sig_t.^-0.5; % Alon 23.6.08
% morlet  = normalize'*ones(1,Ntw).*exp(2*i*pi*frequencies'*t).*exp(power((1./sig_t),2)'*-power(t,2));  % Alon 23.6.08
morlet  = exp( 2*1i*pi*frequencies'*t ) .* exp( - sig_t'.^-2 * t.^2 /2);  % Alon 23.6.08
% morlet = single(morlet);
% signal = single(signal);
% Normalize:
for freq = 1:Nf
    morlet(freq,:) = morlet(freq,:) / norm(morlet(freq,wlim(freq,1):wlim(freq,2)));
end

if Normalize
    for freq = 1:Nf
        morlet(freq,:) = morlet(freq,:) / sum(abs(morlet(freq,:))) * 2;
    end
end

% twin = floor(Ntw/2) + (1:Nts);
% calculating coeffs:
wvlt_amp = single(zeros(Nf,Nts, trials));
wvlt_phase = single(zeros(Nf,Nts, trials));

if ShowWB % 17.10.10
    h = waitbar(0,'Calculating wavelets...');   %   1.1.09
end
% 1-D convolution of each single frequency with each single channel:

    
if UseGpu
    morlet = gpuArray(morlet);
    signal = gpuArray(signal);
end

for trial = (1:trials)
    if ShowText
        disp(['Calculating wavelet ' num2str(trial) ' of ' num2str(trials) '.']);
    end
%     disp(['trial ' num2str(trial) ' out of ' num2str(trials)]);
    if ShowWB % 17.10.10
        if trials > 1 % 14.6.11
            waitbar(trial/trials); %   1.1.09
        end
    end

    for freq = (1:Nf)
        if ShowWB % 14.6.11
            if trials == 1
                waitbar(freq/Nf);
            end
        end
        
        if UseGpu
            Cr = conv(real(morlet(freq,wlim(freq,1):wlim(freq,2))),signal(trial,:));
            Ci = conv(imag(morlet(freq,wlim(freq,1):wlim(freq,2))),signal(trial,:));
            a = gather(Cr(win(freq,:)));
            b = gather(Ci(win(freq,:)));
            wvlt_amp(freq,:,trial) = sqrt(a.^2+b.^2);
            wvlt_phase(freq,:,trial) = atan2(b,a);
        else
            wvlt_trial = conv(morlet(freq,wlim(freq,1):wlim(freq,2)),signal(trial,:));
            wvlt_amp(freq,:,trial) = abs(wvlt_trial(win(freq,:)));   % Should be /sqrt(fs)
            wvlt_phase(freq,:,trial) = angle(wvlt_trial(win(freq,:)));
        end
    end
    
end

if ShowWB % 17.10.10
    close(h)    %   1.1.09
end

% Baseline Correction
if bl_range
    range = bl_range(1):bl_range(2);
    wvlt_amp = bsxfun(@minus,wvlt_amp,mean(wvlt_amp(:,range,:),2));
end

end