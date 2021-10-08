function [total_peak_freq,total_peak_values] = windowpeaks(x,M,H,K,fs,npeaks,dft_peak_freq)

% Returns 2 matrices, having the peak frequency and its magnitude for each
% window

    x = x(:);
    N = length(x);   
       
    kuniq = ceil((K+1)/2); 
    L = 1 + floor((N-M)/H); 
    
    f = (0:kuniq-1)*fs/K;
    
    total_peak_freq = zeros(L,npeaks); %preallocating the matrix
    total_peak_values = zeros(L,npeaks);
    
    actidx = zeros(1,npeaks+1);
   
    for ii = 1:npeaks
        [~,idx]=min(abs(f-dft_peak_freq(ii)));
        actidx(ii+1) = idx;
    end
    
    for l=0:L-1        
        xsw = x(1+l*H : M+l*H);
        xxx = fft(xsw,K);
        xxx = abs(xxx(1:kuniq)).^2/kuniq;
        pi = zeros(1,npeaks); %preallocating the vector
        pvi = zeros(1,npeaks);
        
        for ii = 1:npeaks         
            xxi = xxx(actidx(ii)+1:actidx(ii+1)+2);
             [~, indp] = findpeaks(xxi,'NPeaks',1,'SortStr','descend');  
             if isempty(indp)
                 pi(ii) = NaN;
                 pvi(ii) = NaN;
             else
                pi(ii) = f(indp + actidx(ii)+1);  
                pvi(ii) = xxx(f == pi(ii));
             end
        end
        total_peak_freq(l+1,:) = pi;
        total_peak_values(l+1,:) = pvi;     
    end 
end

