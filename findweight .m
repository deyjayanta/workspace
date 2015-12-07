function y=findweight(nominal,F,m,harmonics)
%  m=10*log10(P);
 for i=1:length(harmonics)
    sig=sum(sum(m(harmonics(i)-i:harmonics(i)+i,:)));
    noise1=sum(sum(m(harmonics(i)-2*i:harmonics(i)-i,:)));
    noise2=sum(sum(m(harmonics(i)+i:harmonics(i)+2*i,:)));
    noise=noise1+noise2;
    snr(i)=sig/noise;
 end
 snr=snr/max(snr);
 weight=snr/sum(snr);
 y=weight;
end