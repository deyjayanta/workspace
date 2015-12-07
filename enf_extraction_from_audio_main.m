
%close all
%clear all
 clc;
 clear all;
 patha={'F:\Level-4 Term-2\IEEE SP Cup\Data\Train_Grid_A_A1.wav'...
        'F:\Level-4 Term-2\IEEE SP Cup\Data\Train_Grid_A_A2.wav'...
        ...
        'F:\Level-4 Term-2\IEEE SP Cup\Data\Train_Grid_B_A1.wav'...
        'F:\Level-4 Term-2\IEEE SP Cup\Data\Train_Grid_B_A2.wav'...
        ...
        'F:\Level-4 Term-2\IEEE SP Cup\Data\Train_Grid_C_A1.wav'...
        'F:\Level-4 Term-2\IEEE SP Cup\Data\Train_Grid_C_A2.wav'...
        ...
        'F:\Level-4 Term-2\IEEE SP Cup\Data\Train_Grid_D_A1.wav'...
        'F:\Level-4 Term-2\IEEE SP Cup\Data\Train_Grid_D_A2.wav'...
        ...
        'F:\Level-4 Term-2\IEEE SP Cup\Data\Train_Grid_E_A1.wav'...
        'F:\Level-4 Term-2\IEEE SP Cup\Data\Train_Grid_E_A2.wav'...
        ...
         'F:\Level-4 Term-2\IEEE SP Cup\Data\Train_Grid_F_A1.wav'...
        'F:\Level-4 Term-2\IEEE SP Cup\Data\Train_Grid_F_A2.wav'...
        ...
        'F:\Level-4 Term-2\IEEE SP Cup\Data\Train_Grid_G_A1.wav'...
        'F:\Level-4 Term-2\IEEE SP Cup\Data\Train_Grid_G_A2.wav'...
        ...
        'F:\Level-4 Term-2\IEEE SP Cup\Data\Train_Grid_H_A1.wav'...
        'F:\Level-4 Term-2\IEEE SP Cup\Data\Train_Grid_H_A2.wav'...
        ...
        'F:\Level-4 Term-2\IEEE SP Cup\Data\Train_Grid_I_A1.wav'...
        'F:\Level-4 Term-2\IEEE SP Cup\Data\Train_Grid_I_A2.wav'...
        };
%  close all;
for j=1:length(patha)
% Y= audioread('D:\sp_cup\Grid_A\Audio_recordings\Train_Grid_A_A1.wav');
% s=datasets( 'audio',j )
% s='Train_Grid_D_A2.wav';
% info=audioinfo(s);
Y= wavread(char(patha(j)));
% info=audioinfo('D:\sp_cup\Grid_B\Power_recordings\Train_Grid_B_P4.wav');
% determining the nominal frequency of grid 
Fs = 1000;

% Y= audioread(s);
% Y= audioread('D:\sp_cup\Grid_B\Power_recordings\Train_Grid_B_P4.wav');
if length(Y)<1000000
    x=Y;
else
x=Y(1:1000000);
end
% x=Y;
[STFT,F,T,P] = spectrogram(x,250,200,1000,Fs);
m=10*log10(P);
%  surf(T,F,m,'edgecolor','none'); axis tight;
% view(0,90);
a=sum(m(50,:))+sum(m(100,:))+sum(m(150,:))+sum(m(200,:))+sum(m(250,:))+sum(m(300,:))+sum(m(350,:));
b=sum(m(60,:))+sum(m(120,:))+sum(m(180,:))+sum(m(240,:))+sum(m(300,:))+sum(m(360,:))+sum(m(420,:));
if a>b
    c=50
    z=sum(m(35,:))+sum(m(70,:))+sum(m(130,:))+sum(m(180,:))+sum(m(230,:))+sum(m(320,:))+sum(m(370,:));
    dec=a/z;
else 
    c=60
    z=sum(m(40,:))+sum(m(80,:))+sum(m(140,:))+sum(m(200,:))+sum(m(260,:))+sum(m(320,:))+sum(m(390,:));
    dec=b/z;
end
% abs(a-b)
if dec>=.85
    sig='audio'
else
    sig='power'
end
% Y1=emd(Y);Y2=(Y1(2,:));
nfft=32768;
% nfft=36000;
harmonic_multiples=1:(ceil(Fs/(2*c))-1);
% harmonic_multiples=1:1;

duration=1;
frame_size_secs=5;%100 samples window(50 hz)
overlap_amount_secs=3;%60 samples overlapping

strip_index=1;
tol=15;
nominal=c;


width_band=[1 3 10];%% 1 %% H,B=3

% for i=3:3
% if i==1
%     width_signal=0.02;
% elseif i==2
%     width_signal=0.06;
% elseif i==3
    width_signal=[0.1 1 8];%%0.1 H,B=1
% elseif i==4
%     width_signal=0.5;
% elseif i==5
%     width_signal=1;        
% end

[enf,type(j)]=findenfhampel(Y, Fs, harmonic_multiples,strip_index, duration, frame_size_secs,overlap_amount_secs, nfft, nominal, width_signal, width_band,tol );
ENF_audio(j,:)=enf';

% figure();
% plot(enf);
% figure
% if i==1
%     plot(enf);
%     hold on;
% elseif i==2
%     plot(enf,'r');
%     hold on;
% elseif i==3
%     plot(enf,'g');
%     hold on;
% elseif i==4
%     plot(enf,'c');
%     hold on;
% elseif i==5
%     plot(enf,'k');
%     hold on;       
% end

end
type;

% figure(1);plot(enfb)
% b=enfb-mean(enfb);hold on;
% figure(2);plot(abs(fft(b)))
% e2=enfb;hold on;

%  R = corrcoef(e1,e2)
% R = corrcoef(enfb(1,:),enfb(2,:))
% R = corrcoef(enfb(1,:),enfb(3,:))
% R = corrcoef(enfb(1,:),enfb(4,:))
% % R = corrcoef(enfb(1,:),enfb(5,:))
% R = corrcoef(enfb(2,:),enfb(3,:))
% R = corrcoef(enfb(2,:),enfb(4,:))
% % R = corrcoef(enfb(2,:),enfb(5,:))
% R = corrcoef(enfb(3,:),enfb(4,:))
% R = corrcoef(enfb(3,:),enfb(5,:))
% R = corrcoef(enfb(4,:),enfb(5,:))