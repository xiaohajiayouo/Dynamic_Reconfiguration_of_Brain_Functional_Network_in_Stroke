function FCM = sliding_window_FC(data,window,step,val )

[nobs, nvar] = size(data);

if window == nobs
    slides=1;
else
    slides=floor((nobs-window)/step)+1; % how many windows
end

FCM = zeros(slides,nvar,nvar);

num0 = ceil(log10(slides))+2;
FC = cell(slides,1);  
t1=1-step;
t2=window-step;
%sliding window 

for k=1:slides
    t1=t1+step;
    t2=t2+step;
    disp([t1 t2]);
    dat = data(t1:t2,:);
    A = corr(dat);
    M = atanh(A - diag(diag(A)));
    if val == 'pos'
        M(M<0) = 0;
    end
    if val == 'abs'
        M = abs(M);
    end
    FCM(k,:,:) = M;
end

end