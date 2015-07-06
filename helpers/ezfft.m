function [x, Y]=ezfft(X,varargin)
% [xf YF]=ezfft(X,['Fs','xf','krn','method']


Fs = 1;  % Sampling frequency
wndw=[];
method='fft';
ax=[];
do_plot=false;

overridedefaults(who,varargin);


num_samp=size(X,2);
num_trials=size(X,1);

NFFT = 2^nextpow2(num_samp); % Next power of 2 from length of y
if isempty(wndw)
    wndw=[0 NFFT/2+1];
end


switch method
    case 'fft'
        Y=nan(num_trials,NFFT);
        for rx=1:num_trials
            Y(rx,:) = abs(fft(X(rx,:),NFFT))/num_samp;
            
        end
        x = Fs/2*linspace(0,1,NFFT/2+1);
    case 'welch'
        Y=nan(num_trials, 2^(ceil(log2(num_samp/9)))+1);
        for rx=1:num_trials
            [t, W]=pwelch(X(rx,:),[],[],[],Fs);
            Y(rx,:)=t';
        end
        x=W';
        
        
end


Y=Y(:,(x>=wndw(1) & x<=wndw(2)));
x=x(x>=wndw(1) & x<=wndw(2));
% Plot single-sided amplitude spectrum.
if do_plot
    if isempty(ax)
        clf; ax=axes;
    else
        set(ax,'NextPlot','add');
    end
    
    plot(ax,x,log10(mean(Y,1)))
    title('Single-Sided Amplitude Spectrum of y(t)')
    xlabel('Frequency (Hz)')
    ylabel('log_10|Y(f)|')
end