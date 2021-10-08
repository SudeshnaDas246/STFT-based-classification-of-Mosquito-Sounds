%clear, clc, close all
close all
format short g

% dir = ["\A\A1","\A\A2","\B\B1","\C\C1","\C\C2","\D\D1","\D\D2","\E\E1","\E\E2"];
% for ii = 1:length(dir)
[x, fs] = audioread(".\Mosquito sound files_13.09.2021\5\5M.wav");

x = x(:, 1);                        % get the first channel

wlen = 1024;                        % window length (recomended to be power of 2)
hop = wlen/4;                       % hop size (recomended to be power of 2)
nfft = 1024;                        % number of fft points (recomended to be power of 2)

win = blackman(wlen, 'periodic');

[S, f, t] = trystft(x, win, hop, nfft, fs);

C = sum(win)/wlen;

S = abs(S)/wlen/C;

% if rem(nfft, 2)                     % odd nfft excludes Nyquist point
%     S(2:end, :) = S(2:end, :).*2;
% else                                % even nfft includes Nyquist point
%     S(2:end-1, :) = S(2:end-1, :).*2;
% end

S = 20*log10(S + 1e-6);

figure(1);
surf(t,f,S)
shading interp
axis tight
view(0, 90)
xlabel('Time, s')
ylabel('Frequency, Hz')
title('Amplitude spectrogram of the signal')
axis([0 inf 0 4000])
hcol = colorbar;
ylabel(hcol, 'Magnitude, dB')

y = fft(x);
n = length(y);          % number of samples
ff = (0:n-1)*(fs/n);     % frequency range
power = abs(y).^2/n;    % power of the DFT

ind = round(4000/(fs/n));
ff = ff(1:ind);
power = power(1:ind);
                        
sgf_sm = sgolayfilt(power, 1, 3501);                        

runavg = movmean(sgf_sm,2000);
% runavg = sgf_sm;
[peak_values,indexes_of_peaks] = findpeaks(runavg,'Npeaks',5,'MinPeakDistance',5000,'SortStr','descend');

% if length(peak_values) == 5
%     [peak_values,indexes_of_peaks] = checkpeaks (runavg,peak_values,indexes_of_peaks);
% end
figure(2);plot(ff,runavg,'-r',ff(indexes_of_peaks),peak_values,'b.', 'MarkerSize', 20,'LineWidth',1);
xlabel('Frequency, Hz')
ylabel('Power')
title('Power spectrum of the signal')


sorted_index_of_peaks = sort(indexes_of_peaks);
indexes_of_valleys = [1];

for jj = 1:length(sorted_index_of_peaks)-1
   indexes_of_valleys = [indexes_of_valleys, find(runavg == min(runavg(sorted_index_of_peaks(jj):sorted_index_of_peaks(jj+1))))];
end

try
    indexes_of_valleys = [indexes_of_valleys, find(runavg == min(runavg(sorted_index_of_peaks(end):sorted_index_of_peaks(end)+10000)))];
catch
    warning('Lack of 5 distinct peaks');
    indexes_of_valleys = [indexes_of_valleys, find(runavg == min(runavg(sorted_index_of_peaks(end):length(runavg))))];    
end

% % figure(2);hold on;plot(ff(indexes_of_valleys),runavg(indexes_of_valleys),'b.', 'MarkerSize', 20);title("peaks")

fwhm_indexes = [];
fwhm_start_freq = [];
fwhm_end_freq = [];
fwhm = [];
for jj = 1:length(sorted_index_of_peaks)
    mini = max(runavg(indexes_of_valleys(jj)),runavg(indexes_of_valleys(jj+1)));
    halfmax = (runavg(sorted_index_of_peaks(jj))+mini)/2;
    index1 = find( runavg(indexes_of_valleys(jj):sorted_index_of_peaks(jj))  >= halfmax, 1, 'first') + indexes_of_valleys(jj);
    index2 = find( runavg(sorted_index_of_peaks(jj):indexes_of_valleys(jj+1))  <= halfmax, 1, 'first') + sorted_index_of_peaks(jj);
    fwhm_indexes = [fwhm_indexes,index1,index2];
    fwhm_start_freq = [fwhm_start_freq,ff(index1)];
    fwhm_end_freq = [fwhm_end_freq,ff(index2)];
    fwhm = [fwhm,ff(index2)-ff(index1)];
end

figure(2);hold on;plot(ff(fwhm_indexes),runavg(fwhm_indexes),'k.','MarkerSize', 20);


sorted_peak_freq = sort(ff(indexes_of_peaks));
sorted_peak_values = zeros(1,length(sorted_peak_freq));
for jj = 1:length(sorted_peak_freq)
    sorted_peak_values(jj) = runavg(ff == sorted_peak_freq(jj));
end

npeaks = 5;

[window_peak_freqs,window_peak_values] = windowpeaks(x, length(win), hop, nfft, fs, npeaks,sorted_peak_freq);

mean_window_peak_freqs = mean(window_peak_freqs,'omitnan');
mean_window_peak_values = mean(window_peak_values,'omitnan');

fft_peak_values = zeros(1,5);
for kk = 1:length(sorted_peak_freq)
    fft_peak_values(kk) = sgf_sm(ff == sorted_peak_freq(kk));
end

% convert to dB
sorted_peak_values = 20*log10(sorted_peak_values + 1e-6);
mean_window_peak_values = 20*log10(mean_window_peak_values + 1e-6);
fft_peak_values = 20*log10(fft_peak_values + 1e-6);

% for kk = 1:length(sorted_peak_freq)
%     yp = sorted_peak_freq(kk);
%     figure(1);hold on;
%     yline(yp,'--r','LineWidth',1);
% end

for kk = 1:length(fwhm_indexes)
    yp = ff(fwhm_indexes(kk));
    figure(1);hold on;
    yline(yp,'--k','LineWidth',1);
end

arr = [sorted_peak_freq',mean_window_peak_freqs',fwhm',fwhm_start_freq',fwhm_end_freq',sorted_peak_values',mean_window_peak_values'];

t = array2table(arr);
t.Properties.VariableNames(1:7) = {'FFT peak freq','Window peak freq','Full Width Half Maximum','FWHM starting freq','FWHM ending freq','FFT peak values', 'Window peak values'};
disp(t);
% % filename = "testfile.csv";
% % f = fopen(filename,'w');
% % fclose(f);
% % writetable(t,filename);
% 
% % end