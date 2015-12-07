function F = power_enf(filename,data_win)
    
    Y = audioread(filename);
    lngth = length(Y);
    data_win = data_win*1000;
    F=[];
    for i=1:((lngth/data_win)-1)
       
        p = Y(data_win*(i-1)+1:data_win*(i-1)+data_win);
        
        p1=fft(p',2^nextpow2(data_win));
        f = (0:length(p1)-1)*1000/length(p1);
        
        w = f>=46 & f<=64;
        p1 = p1.*w;
        
        [~,m]=max(log(abs(p1).^2));
        alpha=20*log(abs(p1(m-1)));
        beta=20*log(abs(p1(m)));
        lambda=20*log(abs(p1(m+1)));
        m1=.5*(alpha-lambda)*(f(m+1)-f(m))/(alpha-2*beta+lambda); % for quadratic interpolation
        f(m)=f(m)+m1;
        F=[F f(m)];
    end
end