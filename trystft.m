function [X,f,t] = trystft(x,w,H,K,fs)

% Returns the STFT matrix and the corresponding frequency and time vectors

    x = x(:);
    N = length(x);   
    M = length(w);
    
    kuniq = ceil((K+1)/2); %The number of unique fft points
    L = 1 + floor((N-M)/H); %The number of the successive signal frames across the length of the analyzed time sequence
    
    X = zeros(kuniq,L);
    
    for l=0:L-1
        ssw = x(1+l*H : M+l*H).*w;
        xl = fft(ssw,K);
        X(1:kuniq , 1+l) = xl(1:kuniq);
    end 
    
    t = (M/2 : H : M/2+(L-1)*H)/fs;
    f = (0:kuniq-1)*fs/K;
    
end

