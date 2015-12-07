function [enf,sig,c] = enf_extract_from_audio_and_power(file,win,overlap)

Fs = 1000;

Y= audioread(file);

if length(Y)<1000000
    x=Y;
else
x=Y(1:1000000);
end
% x=Y;
[STFT,F,T,P] = spectrogram(x,250,200,1000,Fs);
 m=10*log10(P);
         %surf(T,F,m,'edgecolor','none'); axis tight;
        %view(0,90);
        a=sum(m(50,:))+sum(m(100,:))+sum(m(150,:))+sum(m(200,:))+sum(m(250,:))+sum(m(300,:))+sum(m(350,:));
        a1=sum(P(50,:))+sum(P(100,:))+sum(P(150,:))+sum(P(200,:))+sum(P(250,:))+sum(P(300,:))+sum(P(350,:));
        b1=sum(P(60,:))+sum(P(120,:))+sum(P(180,:))+sum(P(240,:))+sum(P(300,:))+sum(P(360,:))+sum(P(420,:));
       
        b=sum(m(60,:))+sum(m(120,:))+sum(m(180,:))+sum(m(240,:))+sum(m(300,:))+sum(m(360,:))+sum(m(420,:));
        if a>b
            c=50;
            x1=sum(sum(P));
            dec=x1/a1;
      else 
            c=60;
x1=sum(sum(P));
dec=x1/b1;

        end
 if dec>8.5
       sig='audio';
 else
       sig='power';
 end



nfft=32768;

harmonic_multiples=1:(ceil(Fs/(2*c))-1);
% harmonic_multiples=1:1;


duration = 5;
frame_size_secs=win;%100 samples window(50 hz)
overlap_amount_secs=overlap;%60 samples overlapping



strip_index=1;
tol=15;
nominal=c;
width_band=[1 3 8];


    width_signal=[0.1 1 5];

if sig == 'audio'
   
   [enf,~]=findenfhampel(Y, Fs, harmonic_multiples,strip_index, duration, frame_size_secs,overlap_amount_secs, nfft, nominal, width_signal, width_band,tol );
else
    enf = power_enf(file,2);
end


end