function [ ENF,type ] = findenfhampel( signal, fs, harmonic_multiples,strip_index, duration, frame_size_secs, overlap_amount_secs, nfft, nominal, width_signal, width_band,tol)

for i=1:length(width_signal)
harmonics=nominal*harmonic_multiples;
[ weights] = computeCombiningWeights( signal, fs, harmonics, duration, frame_size_secs, overlap_amount_secs, nfft, nominal, width_signal(i), width_band(i) );

[ strips, frequency_support ] = computeSpectrogramStrips( signal, frame_size_secs, overlap_amount_secs, nfft, fs, nominal, harmonic_multiples, width_band(i) );
% [ pseudo_spectra_strips, pseudo_frequency_support ] = computePseudoSpectrumStrips( signal, filters, frame_size_secs, overlap_amount_secs, nfft, fs, nominal, harmonic_multiples, width_band (i) );

[ OurStrip_Cell, initial_frequency ] = computeCombinedSpectrum( strips,weights, duration, frame_size_secs,overlap_amount_secs, strip_index, frequency_support );

[ENF] = computeENFfromCombinedStrip( OurStrip_Cell, initial_frequency, fs, nfft);

% [ OurStrip_Cell, initial_frequency ] = computeCombinedSpectrum( pseudo_spectra_strips,weights, duration, frame_size_secs, strip_index, pseudo_frequency_support );
% [pseudo_ENF] = computeENFfromCombinedStrip( OurStrip_Cell, initial_frequency, fs, nfft);

% figure()
% plot(ENF);
% hold on;


ENF=hampelfilt(1:length(ENF),ENF);

ENF=smooth(ENF);
E(i,:)=ENF;
% 
% plot(ENF,'r');

end

D=mean(abs(diff(E,1,2)),2)
ENF=E(min(find(D==min(D))),:);
type=min(find(D==min(D)));
ENF(find(abs(ENF-nominal)>tol))=nominal;

% SD=std(E,0,2);
% ENF=E(min(find(SD==min(SD))),:);
% type=min(find(SD==min(SD)))

% MD=mad(E,0,2);
% ENF=E(min(find(MD==min(MD))));
% type=min(find(MD==min(MD)))


% figure(1)
% plot(ENF,'r');
% hold on;
% plot(pseudo_ENF,'k');
% 
% pseudo_ENF(find(abs(pseudo_ENF-nominal)>tol))=nominal;
% pseudo_ENF=hampelfilt(1:length(pseudo_ENF),pseudo_ENF);
% pseudo_ENF=smooth(pseudo_ENF);
% 
% plot(pseudo_ENF,'c');
% 
% ENF=(ENF+pseudo_ENF)/2;
% 
% plot(ENF,'o');


% ENF=medfilt1(ENF,3);
% d=abs(diff(ENF));
% di=ceil(d/min(d));
% md=mode(di)*min(d);
% mn=mean(ENF);
% for i=1:length(ENF)
%     if i>2 && i<length(ENF)-1
%         t=mean([abs(ENF(i)-ENF(i-2)) abs(ENF(i)-ENF(i-1)) abs(ENF(i)-ENF(i+1)) abs(ENF(i)-ENF(i+2))]);
%         if t/abs(ENF(i+1)-ENF(i-1))>tol          
%             ENF(i)=ENF(i)*.02+ENF(i-1)*0.49+ENF(i+1)*0.49;
%         end
%     elseif i<3
%         t=abs(ENF(i)-mn);
%         if t/md >tol
%            ENF(i)=ENF(i)*.03+ENF(i+1)*0.18+ENF(i+2)*0.19+ENF(i+3)*0.20+ENF(i+4)*0.20+ENF(i+4)*0.20;
%         end
%     else
%         t=mean([abs(ENF(i)-ENF(i-2)) abs(ENF(i)-ENF(i-1))]);
%         if t/abs(ENF(i-1)-ENF(i-2)) >tol
%            ENF(i)=ENF(i)*.02+ENF(i-2)*0.49+ENF(i-1)*0.49;
%         end
%     end   
% 
%     
% end
end

%% Spectrum Combining for ENF Signal Estimation
% Adi Hajj-Ahmad, Student Member, IEEE,RaviGarg, Student Member, IEEE,and MinWu, Fellow, IEEE

function [ weights] = computeCombiningWeights( signal0, fs, harmonics, duration, frame_size_secs, overlap_amount_secs, nfft, nominal, width_signal, width_band )
%COMPUTECOMBININGWEIGHTS Summary of this function goes here
%   This function computes the Combining Weights for the Spectrum Combining
%   approach for ENF signal estimation
%   Takes as input:
%   -> signal0: ENF-containing signal.
%   -> fs: sampling frequency of 'signal0'.
%   -> harmonics: harmonics for which we want to compute weights around.
%   -> duration: time duration for which we compute a certain weight, e.g.
%   30 min
%   -> frame_size_secs: size of time-frame for which one instantaneous ENF
%   point is estimated, e.g. 5 sec.
%   -> overlap_amount_secs: nb of seconds of overlap between time-frames
%   -> nfft: Nb of points for FFT computation, e.g. 32768 for frequency
%   resolution of ~0.03Hz when fs = 1000Hz.
%   -> nominal: nominal frequency, e.g. 60Hz for US grids.
%   -> width_signal: half the width of the band which we consider contains
%   the ENF fluctutations around the nominal value, e.g. 0.02Hz for US ENF.
%   -> width_band: half the width of the band abour the nominal value,
%   which we are considering for SNR computations, e.g. 1Hz for US ENF.
%   Gives output:
%   -> spectro_strips: a Matlab Cell of size equal to the number of
%   harmonic multiples, each component contains a spectrogram strip
%   centered at one of the harmonic multiples.
%   Gives output:
%   -> weights: combining weights computed.

% setting up the variables
nb_durations = ceil(length(signal0)/(duration*60*fs));
frame_size = floor(fs*frame_size_secs);
overlap_amount = floor(fs*overlap_amount_secs);
shift_amount = frame_size - overlap_amount;
nb_of_harmonics = length(harmonics);
harmonic_multiples = harmonics/nominal;
starting_freq = nominal - width_band;
center_freq = nominal;
init_first_value = nominal - width_signal;
init_second_value = nominal + width_signal;
weights = zeros(nb_of_harmonics,nb_durations);
inside_mean = zeros(nb_of_harmonics,nb_durations);
outside_mean = zeros(nb_of_harmonics,nb_durations);
total_nb_frames = 0;
All_Strips_Cell =  cell(nb_durations, 1);

for dur = 1:nb_durations
    
    x = signal0( (dur -1)*duration*60*fs + 1: min(end, dur*duration*60*fs + overlap_amount));
    
    %% getting the spectrogram 
    nb_of_frames = ceil((length(x) - frame_size + 1)/shift_amount);
    total_nb_frames = total_nb_frames + nb_of_frames;
    P = zeros(nfft/2 + 1, nb_of_frames);
    starting = 1;
    for frame = 1:nb_of_frames
        ending = starting + frame_size - 1;
        signal = x(starting:ending);
        [S, F, T ,P(:, frame)] = spectrogram(signal, frame_size, 0, nfft, fs);
        starting = starting + shift_amount;
    end
    
    %% getting the harmonic strips
    
    width_init = findClosest(F, center_freq) - findClosest(F, starting_freq);
    HarmonicStrips = zeros(width_init*2*sum(harmonic_multiples) , size(P, 2));
    FreqAxis = zeros(width_init*2*sum(harmonic_multiples), 1);
    resolution = F(2) - F(1);

    starting = 1;
    starting_indices = zeros(nb_of_harmonics,1);
    ending_indices = zeros(nb_of_harmonics,1);
    for k = 1:nb_of_harmonics
        starting_indices(k) = starting;
        width = width_init*harmonic_multiples(k);
        ending = starting + 2*width -1;
        ending_indices(k) = ending;
        tempFreqIndex = round(harmonics(k)/resolution);
        
        %%
        sss=size(P);
        if tempFreqIndex + width>sss(1)
            ttt=sss(1)-(ending-starting);
            st=sss(1);
        else
            ttt=tempFreqIndex - width + 1;
            st=tempFreqIndex + width;
        end
        HarmonicStrips(starting:ending, :) = P(ttt:st, :);
        FreqAxis(starting:ending) = F(ttt:st, :);
%         HarmonicStrips(starting:ending, :) = P((tempFreqIndex - width + 1):(tempFreqIndex + width), :);
 %       FreqAxis(starting:ending) = F((tempFreqIndex - width + 1):(tempFreqIndex + width));
        starting = ending + 1;
    end
    
    All_Strips_Cell{dur} = HarmonicStrips;
    %% getting the weights
    
    for k = 1:nb_of_harmonics
        currStrip = HarmonicStrips(starting_indices(k):ending_indices(k),:);
        freq_axis = FreqAxis(starting_indices(k):ending_indices(k));
        first_value = init_first_value*harmonic_multiples(k);
        second_value = init_second_value*harmonic_multiples(k);
        first_index = findClosest(freq_axis, first_value);
        second_index = findClosest(freq_axis, second_value);
        inside_strip = currStrip(first_index:second_index, :);
        inside_mean(k, dur) = mean(mean(inside_strip));
        outside_strip = currStrip([1:first_index-1, second_index +1:end], :);
        outside_mean(k, dur) = mean(mean(outside_strip));
        if inside_mean(k, dur) < outside_mean(k, dur)
            weights(k, dur) = 0;
        else
            weights(k, dur) = inside_mean(k, dur)/outside_mean(k, dur);
        end
        
    end
    
    %% normalizing the weights
    sum_weights = sum(weights(:, dur));
    for k = 1:nb_of_harmonics
        weights(k, dur) = 100*weights(k, dur)/sum_weights;
    end
end

end


function [ spectro_strips, frequency_support ] = computeSpectrogramStrips( signal, frame_size_secs, overlap_amount_secs, nfft, fs, nominal, harmonic_multiples, width_band )
%COMPUTESPECTROGRAMSTRIPS Summary of this function goes here
%   This function generates the spectrogram strips needed for ENF signal estimation
%   Takes as input: 
%   -> signal: ENF-containing signal.
%   -> frame_size_secs: size of time-frame for which one instantaneous ENF
%   point is estimated, e.g. 5 sec.
%   -> overlap_amount_secs: nb of seconds of overlap between time-frames
%   -> nfft: Nb of points for FFT computation, e.g. 32768 for frequency
%   resolution of ~0.03Hz when fs = 1000Hz.
%   -> fs: sampling frequency of 'signal'.
%   -> nominal: nominal frequency of ENF signal, e.g. 60Hz in US.
%   -> harmonic_multiples: harmonic multiples for which we want to compute
%   the spectrogram strips for, e.g. 1:8 for 60, 120, 180, ... , 480Hz.
%   -> width_band: half the desired width of the strips, about the nominal
%   frequency.
%   Gives output:
%   -> spectro_strips: a Matlab Cell of size equal to the number of
%   harmonic multiples, each component contains a spectrogram strip
%   centered at one of the harmonic multiples.
%   -> frequency_support: index of starting and ending frequencies of each
%   strip.

% setting up the variables
nb_harmonics = length(harmonic_multiples);
spectro_strips = cell(nb_harmonics, 1);
frame_size = frame_size_secs*fs;
overlap_amount = overlap_amount_secs*fs;
shift_amount = frame_size - overlap_amount;
len_sig = length(signal);
nb_frames = ceil((len_sig - frame_size + 1)/shift_amount);

% collecting the full spectrogram in P, and the full frequency axis in F
starting = 1;
P = zeros(nfft/2 + 1, nb_frames);
for frame = 1:nb_frames
    ending = starting + frame_size - 1;
    x = signal(starting:ending);
   
    [~, F, ~, P(:, frame)] = spectrogram(x, frame_size, 0, nfft, fs);
    
    starting = starting + shift_amount;
end


% choosing the strips that we need, and setting up 'frequency_support'.
first_index = findClosest(F, nominal - width_band);
second_index = findClosest(F, nominal + width_band);
frequency_support = zeros(nb_harmonics, 2);
for k = 1:nb_harmonics
    starting = first_index*harmonic_multiples(k);
    ending = second_index*harmonic_multiples(k);
%     k
%     size(P)
%     starting
%     ending
    sss=size(P);
    if ending>sss(1)
        temp=ending-starting;
        ending=sss(1);
        starting=ending-temp;
    end
        
    
    spectro_strips{k} = P(starting:ending, :);
    
    frequency_support(k, 1) = F(starting);
    frequency_support(k, 2) = F(ending);
end

end

function [ OurStrip_Cell, initial_frequency ] = computeCombinedSpectrum( strips, weights, duration, frame_size_secs,overlap_amount_secs, strip_index, frequency_support )
%COMPUTECOMBINEDSPECTRUM Summary of this function goes here
%   This function takes in spectrogram strips, or pseudo-spectrum strips
%   and combines them according to the combining weights.
%   Takes as input:
%   -> strips: Cell containing strips around different harmonics.
%   -> weights: combining weights.
%   -> duration: time-duration for which a certain weight is computed.
%   -> frame_size_secs: size of time-frame for which one instantaneous ENF
%   point is estimated, e.g. 5 sec.
%   -> strip_index: index of the strip in 'strips' whose width will be 
%   the width we give to all the strips when we combine them.
%   -> frequency_support: beginning and ending frequencies for each strip
%   Gives as output:
%   -> OurStrip_cell: Cell of size corresponding to the number of
%   durations, containing the combined strip for each duration.
%   -> initial_frequency: starting frequency of combined strip

% setting up the variables
nb_durations = size(weights, 2);
nb_frames = size(strips{1}, 2);
%%
nb_frames_per_dur = floor(duration*60/(frame_size_secs-overlap_amount_secs));
%% nb_frames_per_dur = duration*60/frame_size_secs;
strip_width = size(strips{strip_index}, 1);
OurStrip_Cell = cell(nb_durations, 1);
nb_signals = length(strips);
initial_frequency = frequency_support(strip_index, 1);

% combining the strips, taking each duration at a time, as each duration
% has a different set of weights.
begin = 1;
for dur = 1:nb_durations
    nb_frames_left = nb_frames - (dur-1)*nb_frames_per_dur;
    OurStrip = zeros(strip_width, min(nb_frames_per_dur, nb_frames_left));
    endit = begin + size(OurStrip,2) -1;
    for harm = 1:nb_signals
        tempStrip = strips{harm}(:, begin:endit);
        for frame = 1:size(OurStrip, 2)
            tempo = imresize(tempStrip(:, frame), [strip_width, 1], 'bilinear');
            tempo = 100*tempo/max(tempo);
            OurStrip(:, frame) = OurStrip(:, frame) + weights(harm, dur)*tempo;
        end
    end
    OurStrip_Cell{dur} = OurStrip;
    begin = endit +1 ;
end
end

function [ pseudo_spectra_strips, frequency_support ] = computePseudoSpectrumStrips( signal, filters, frame_size_secs, overlap_amount_secs, nfft, fs, nominal, harmonic_multiples, width_band  )
%COMPUTEPSEUDOSPECTRUMSTRIPS Summary of this function goes here
%   This function generates the pseudo-spectrum strips needed for ENF signal estimation
%   Takes as input: 
%   -> signal: ENF-containing signal.
%   -> filters: Matlab Cell containing band-pass filters centered about the
%   desired harmonics.
%   -> frame_size_secs: size of time-frame for which one instantaneous ENF
%   point is estimated, e.g. 5 sec.
%   -> overlap_amount_secs: nb of seconds of overlap between time-frames
%   -> nfft: Nb of points for FFT computation, e.g. 32768 for frequency
%   resolution of ~0.03Hz  when fs = 1000Hz.
%   -> fs: sampling frequency of 'signal.
%   -> nominal: nominal frequency of ENF signal, e.g. 60Hz in US.
%   -> harmonic_multiples: harmonic multiples for which we want to compute
%   the spectrogram strips for, e.g. 1:8 for 60, 120, 180, ... , 480Hz.
%   -> width_band: half the desired width of the strips, about the nominal
%   frequency.
%   Gives output:
%   -> pseudo_spectra_strips: a Matlab Cell of size equal to the number of
%   harmonic multiples, each component contains a pseudospectrum strip
%   centered at one of the harmonic multiples.
%   -> frequency_support: index of starting and ending frequencies of each
%   strip.

% setting up the variables
nb_harmonics = length(filters);
pseudo_spectra_strips = cell(nb_harmonics, 1);
signals = filter_signals(signal, filters, 1:nb_harmonics);
frame_size = frame_size_secs*fs;
overlap_amount = overlap_amount_secs*fs;
shift_amount = frame_size - overlap_amount;
len_sig = length(signal);
nb_frames = ceil((len_sig - frame_size + 1)/shift_amount);
frequency_support = zeros(nb_harmonics, 2);

% computing full pseudo-spectrum of each filtered signal, then choosing out
% the strips needed.
for k = 1:nb_harmonics
        x = signals{k};
        pseudo_temp = zeros(nfft/2 + 1, nb_frames);
        starting = 1;
        for frame = 1:nb_frames
            ending = starting + frame_size -1;
            [S, f] = pmusic(x(starting:ending), 2, nfft, fs, 50, 49);
            pseudo_temp(:, frame) = S;   
            starting = starting + shift_amount;
        end
        
        first_index = findClosest(f, nominal - width_band);
        second_index = findClosest(f, nominal + width_band);
        
        starting = first_index*harmonic_multiples(k);
        ending = second_index*harmonic_multiples(k);
        pseudo_spectra_strips{k} = pseudo_temp(starting:ending, :);
        frequency_support(k, 1) = f(starting);
        frequency_support(k, 2) = f(ending);
        clear pseudo_temp
end
end

function [ENF] = computeENFfromCombinedStrip( OurStrip_Cell, initial_frequency, fs, nfft  )
%COMPUTEENFFROMCOMBINEDSTRIP Summary of this function goes here
%   This function computes the ENF signal estimate from the combined strip
%   Takes as input:
%   -> OurStrip_Cell: Cell containing the combined strip for different
%   durations.
%   -> initial_frequency: the frequency to which the first column in the
%   combined strip corresponds to.
%   -> fs: sampling frequency.
%   -> nfft: Nb of points for FFT computation, e.g. 32768 for resolutoon of
%   ~0.03Hz when fs = 1000Hz.

% setting up the variables
nb_durations = length(OurStrip_Cell);
nb_frames_per_dur = size(OurStrip_Cell{1}, 2);
nb_frames = nb_frames_per_dur*(nb_durations - 1) + size(OurStrip_Cell{end},2);
ENF = zeros(nb_frames, 1);

% taking each frame at a time, and using quadratic interpolation to find
% its maximum and thus its dominant ENF frequency

starting = 1;
for dur = 1:nb_durations
    OurStrip_here = OurStrip_Cell{dur};
    nb_frames_here = size(OurStrip_here, 2);
    ending = starting + nb_frames_here - 1;
    ENF_here = zeros(nb_frames_here, 1);
    for frame = 1:nb_frames_here
        power_vector = OurStrip_here(:, frame);
        [~, index] = max(power_vector);
        k_star = QuadInterpFunction(power_vector, index);
        ENF_here(frame) = initial_frequency + fs*k_star/nfft;
    end
    ENF(starting:ending) = ENF_here;
    starting = ending + 1;
end
end

function [ signals ] = filter_signals( signal, filters, indices )
%FILTER_SIGNALS Summary of this function goes here
%   This function filters 'signal' through the filters in the cell
%   structure 'filters' indexed by the values in 'indices'. The output
%   filtered signals are stored in the cell structure 'signals'.

signal = signal(:);
nb_filters_chosen = length(indices);
signals = cell(nb_filters_chosen, 1);
for k = 1:nb_filters_chosen
    h = filters{indices(k)};
    signals{k} = filter(h.Numerator, 1, signal);
end

end

function [ index ] = findClosest( vector, value )
%FINDCLOSEST Summary of this function goes here
%   This function finds the index of closest element in 'vector' to 'value'

index = 1;
for k = 2:length(vector)
    if (abs(vector(k) - value) < abs(vector(k-1) - value))
        index = k;
    else
        break;
    end
end
end

function [ k_star ] = QuadInterpFunction( vector, index )
%QUAD_INTERP_FUNC Summary of this function goes here
%   This function finds the index 'k_star' of the maximum value in 'vector'
%   about 'index', computed using quadratic interpolation

if index == 1
    index = 2;
elseif index == length(vector)
    index = length(vector) - 1;
end
alpha = 20 * log10(abs(vector(index - 1)));
beta = 20* log10(abs(vector(index)));
gamma = 20*log10(abs(vector(index + 1)));
delta = 0.5*(alpha - gamma)/(alpha - 2*beta + gamma);
kmax = index - 1;
k_star = kmax + delta;
end

function [YY, I, Y0, LB, UB, ADX, NO] = hampelfilt(X, Y, DX, T, varargin)
% HAMPEL    Hampel Filter.
%   HAMPEL(X,Y,DX,T,varargin) returns the Hampel filtered values of the 
%   elements in Y. It was developed to detect outliers in a time series, 
%   but it can also be used as an alternative to the standard median 
%   filter.
%
%   References
%   Chapters 1.4.2, 3.2.2 and 4.3.4 in Mining Imperfect Data: Dealing with 
%   Contamination and Incomplete Records by Ronald K. Pearson.
%
%   Acknowledgements
%   I would like to thank Ronald K. Pearson for the introduction to moving
%   window filters. Please visit his blog at:
%   http://exploringdatablog.blogspot.com/2012/01/moving-window-filters-and
%   -pracma.html
%
%   X,Y are row or column vectors with an equal number of elements.
%   The elements in Y should be Gaussian distributed.
%
%   Input DX,T,varargin must not contain NaN values!
%
%   DX,T are optional scalar values.
%   DX is a scalar which defines the half width of the filter window. 
%   It is required that DX > 0 and DX should be dimensionally equivalent to
%   the values in X.
%   T is a scalar which defines the threshold value used in the equation
%   |Y - Y0| > T*S0.
%
%   Standard Parameters for DX and T:
%   DX  = 3*median(X(2:end)-X(1:end-1)); 
%   T   = 3;
%
%   varargin covers addtional optional input. The optional input must be in
%   the form of 'PropertyName', PropertyValue.
%   Supported PropertyNames: 
%   'standard': Use the standard Hampel filter. 
%   'adaptive': Use an experimental adaptive Hampel filter. Explained under
%   Revision 1 details below.
% 
%   Supported PropertyValues: Scalar value which defines the tolerance of
%   the adaptive filter. In the case of standard Hampel filter this value 
%   is ignored.
%
%   Output YY,I,Y0,LB,UB,ADX are column vectors containing Hampel filtered
%   values of Y, a logical index of the replaced values, nominal data,
%   lower and upper bounds on the Hampel filter and the relative half size 
%   of the local window, respectively.
%
%   NO is a scalar that specifies the Number of Outliers detected.
%
%   Examples
%   1. Hampel filter removal of outliers
%       X           = 1:1000;                           % Pseudo Time
%       Y           = 5000 + randn(1000, 1);            % Pseudo Data
%       Outliers    = randi(1000, 10, 1);               % Index of Outliers
%       Y(Outliers) = Y(Outliers) + randi(1000, 10, 1); % Pseudo Outliers
%       [YY,I,Y0,LB,UB] = hampel(X,Y);
%
%       plot(X, Y, 'b.'); hold on;      % Original Data
%       plot(X, YY, 'r');               % Hampel Filtered Data
%       plot(X, Y0, 'b--');             % Nominal Data
%       plot(X, LB, 'r--');             % Lower Bounds on Hampel Filter
%       plot(X, UB, 'r--');             % Upper Bounds on Hampel Filter
%       plot(X(I), Y(I), 'ks');         % Identified Outliers
%
%   2. Adaptive Hampel filter removal of outliers
%       DX          = 1;                                % Window Half size
%       T           = 3;                                % Threshold
%       Threshold   = 0.1;                              % AdaptiveThreshold
%       X           = 1:DX:1000;                        % Pseudo Time
%       Y           = 5000 + randn(1000, 1);            % Pseudo Data
%       Outliers    = randi(1000, 10, 1);               % Index of Outliers
%       Y(Outliers) = Y(Outliers) + randi(1000, 10, 1); % Pseudo Outliers
%       [YY,I,Y0,LB,UB] = hampel(X,Y,DX,T,'Adaptive',Threshold);
%
%       plot(X, Y, 'b.'); hold on;      % Original Data
%       plot(X, YY, 'r');               % Hampel Filtered Data
%       plot(X, Y0, 'b--');             % Nominal Data
%       plot(X, LB, 'r--');             % Lower Bounds on Hampel Filter
%       plot(X, UB, 'r--');             % Upper Bounds on Hampel Filter
%       plot(X(I), Y(I), 'ks');         % Identified Outliers
%
%   3. Median Filter Based on Filter Window
%       DX        = 3;                        % Filter Half Size
%       T         = 0;                        % Threshold
%       X         = 1:1000;                   % Pseudo Time
%       Y         = 5000 + randn(1000, 1);    % Pseudo Data
%       [YY,I,Y0] = hampel(X,Y,DX,T);
%
%       plot(X, Y, 'b.'); hold on;    % Original Data
%       plot(X, Y0, 'r');             % Median Filtered Data
%
%   Version: 1.5
%   Last Update: 09.02.2012
%
%   Copyright (c) 2012:
%   Michael Lindholm Nielsen
%
%   --- Revision 5 --- 09.02.2012
%   (1) Corrected potential error in internal median function.
%   (2) Removed internal "keyboard" command.
%   (3) Optimized internal Gauss filter.
%
%   --- Revision 4 --- 08.02.2012
%   (1) The elements in X and Y are now temporarily sorted for internal
%       computations.
%   (2) Performance optimization.
%   (3) Added Example 3.
%
%   --- Revision 3 --- 06.02.2012
%   (1) If the number of elements (X,Y) are below 2 the output YY will be a
%       copy of Y. No outliers will be detected. No error will be issued.
%
%   --- Revision 2 --- 05.02.2012
%   (1) Changed a calculation in the adaptive Hampel filter. The threshold
%       parameter is now compared to the percentage difference between the
%       j'th and the j-1 value. Also notice the change from Threshold = 1.1
%       to Threshold = 0.1 in example 2 above.
%   (2) Checks if DX,T or varargin contains NaN values.
%   (3) Now capable of ignoring NaN values in X and Y.
%   (4) Added output Y0 - Nominal Data.
%
%   --- Revision 1 --- 28.01.2012
%   (1) Replaced output S (Local Scaled Median Absolute Deviation) with
%       lower (LB) and upper (UB) bounds on the Hampel filter.
%   (2) Added option to use an experimental adaptive Hampel filter.
%       The Principle behind this filter is described below.
%   a) The filter changes the local window size until the change in the 
%       local scaled median absolute deviation is below a threshold value 
%       set by the user. In the above example (2) this parameter is set to 
%       0.1 corresponding to a maximum acceptable change of 10% in the 
%       local scaled median absolute deviation. This process leads to three
%       locally optimized parameters Y0 (Local Nominal Data Reference 
%       value), S0 (Local Scale of Natural Variation), ADX (Local Adapted 
%       Window half size relative to DX).
%   b) The optimized parameters are then smoothed by a Gaussian filter with
%       a standard deviation of DX=2*median(XSort(2:end) - XSort(1:end-1)).
%       This means that local values are weighted highest, but nearby data 
%       (which should be Gaussian distributed) is also used in refining 
%       ADX, Y0, S0.
%   
%   --- Revision 0 --- 26.01.2012
%   (1) Release of first edition.

%% Error Checking
% Check for correct number of input arguments
if nargin < 2
    error('Not enough input arguments.');
end

% Check that the number of elements in X match those of Y.
if ~isequal(numel(X), numel(Y))
    error('Inputs X and Y must have the same number of elements.');
end

% Check that X is either a row or column vector
if size(X, 1) == 1
    X   = X';   % Change to column vector
elseif size(X, 2) == 1
else
    error('Input X must be either a row or column vector.')
end

% Check that Y is either a row or column vector
if size(Y, 1) == 1
    Y   = Y';   % Change to column vector
elseif size(Y, 2) == 1
else
    error('Input Y must be either a row or column vector.')
end

% Sort X
SortX   = sort(X);

% Check that DX is of type scalar
if exist('DX', 'var')
    if ~isscalar(DX)
        error('DX must be a scalar.');
    elseif DX < 0
        error('DX must be larger than zero.');
    end
else
    DX  = 3*median(SortX(2:end) - SortX(1:end-1));
end

% Check that T is of type scalar
if exist('T', 'var')
    if ~isscalar(T)
        error('T must be a scalar.');
    end
else
    T   = 3;
end

% Check optional input
if isempty(varargin)
    Option  = 'standard';
elseif numel(varargin) < 2
    error('Optional input must also contain threshold value.');
else
    % varargin{1}
    if ischar(varargin{1})
        Option      = varargin{1};
    else
        error('PropertyName must be of type char.');
    end
    % varargin{2}
    if isscalar(varargin{2})
        Threshold   = varargin{2};
    else
        error('PropertyValue value must be a scalar.');
    end
end

% Check that DX,T does not contain NaN values
if any(isnan(DX) | isnan(T))
    error('Inputs DX and T must not contain NaN values.');
end

% Check that varargin does not contain NaN values
CheckNaN    = cellfun(@isnan, varargin, 'UniformOutput', 0);
if any(cellfun(@any, CheckNaN))
    error('Optional inputs must not contain NaN values.');
end

% Detect/Ignore NaN values in X and Y
IdxNaN  = isnan(X) | isnan(Y);
X       = X(~IdxNaN);
Y       = Y(~IdxNaN);

%% Calculation
% Preallocation
YY  = Y;
I   = false(size(Y));
S0  = NaN(size(YY));
Y0  = S0;
ADX = repmat(DX, size(Y));

if numel(X) > 1
    switch lower(Option)
        case 'standard'
            for i = 1:numel(Y)
                % Calculate Local Nominal Data Reference value
                % and Local Scale of Natural Variation
                [Y0(i), S0(i)]  = localwindow(X, Y, DX, i);
            end
        case 'adaptive'
            % Preallocate
            Y0Tmp   = S0;
            S0Tmp   = S0;
            DXTmp   = (1:numel(S0))'*DX; % Integer variation of Window Half Size
            
            % Calculate Initial Guess of Optimal Parameters Y0, S0, ADX
            for i = 1:numel(Y)
                % Setup/Reset temporary counter etc.
                j       = 1;
                S0Rel   = inf;
                while S0Rel > Threshold
                    % Calculate Local Nominal Data Reference value
                    % and Local Scale of Natural Variation using DXTmp window
                    [Y0Tmp(j), S0Tmp(j)]    = localwindow(X, Y, DXTmp(j), i);
                    
                    % Calculate percent difference relative to previous value
                    if j > 1
                        S0Rel   = abs((S0Tmp(j-1) - S0Tmp(j))/(S0Tmp(j-1) + S0Tmp(j))/2);
                    end
                    
                    % Iterate counter
                    j   = j + 1;
                end
                Y0(i)   = Y0Tmp(j - 2);     % Local Nominal Data Reference value
                S0(i)   = S0Tmp(j - 2);     % Local Scale of Natural Variation
                ADX(i)  = DXTmp(j - 2)/DX;  % Local Adapted Window size relative to DX
            end
            
            % Gaussian smoothing of relevant parameters
            DX  = 2*median(SortX(2:end) - SortX(1:end-1));
            ADX = smgauss(X, ADX, DX);
            S0  = smgauss(X, S0, DX);
            Y0  = smgauss(X, Y0, DX);
        otherwise
            error('Unknown option ''%s''.', varargin{1});
    end
end

%% Prepare Output
UB      = Y0 + T*S0;            % Save information about local scale
LB      = Y0 - T*S0;            % Save information about local scale
Idx     = abs(Y - Y0) > T*S0;   % Index of possible outlier
YY(Idx) = Y0(Idx);              % Replace outliers with local median value
I(Idx)  = true;                 % Set Outlier detection
NO      = sum(I);               % Output number of detected outliers

% Reinsert NaN values detected at error checking stage
if any(IdxNaN)
    [YY, I, Y0, LB, UB, ADX]    = rescale(IdxNaN, YY, I, Y0, LB, UB, ADX);
end

%% Built-in functions
    function [Y0, S0] = localwindow(X, Y, DX, i)
        % Index relevant to Local Window
        Idx = X(i) - DX <= X & X <= X(i) + DX;

        % Calculate Local Nominal Data Reference Value
        Y0  = median(Y(Idx));
        
        % Calculate Local Scale of Natural Variation
        S0  = 1.4826*median(abs(Y(Idx) - Y0));
    end

    function M = median(YM)
        % Isolate relevant values in Y
        YM  = sort(YM);
        NYM = numel(YM);
        
        % Calculate median
        if mod(NYM,2)   % Uneven
            M   = YM((NYM + 1)/2);
        else            % Even
            M   = (YM(NYM/2)+YM(NYM/2+1))/2;
        end
    end

    function G = smgauss(X, V, DX)
        % Prepare Xj and Xk
        Xj  = repmat(X', numel(X), 1);
        Xk  = repmat(X, 1, numel(X));
        
        % Calculate Gaussian weight
        Wjk = exp(-((Xj - Xk)/(2*DX)).^2);
        
        % Calculate Gaussian Filter
        G   = Wjk*V./sum(Wjk,1)';
    end

    function varargout = rescale(IdxNaN, varargin)
        % Output Rescaled Elements
        varargout    = cell(nargout, 1);
        for k = 1:nargout
            Element     = varargin{k};
            
            if islogical(Element)
                ScaledElement   = false(size(IdxNaN));
            elseif isnumeric(Element)
                ScaledElement   = NaN(size(IdxNaN));
            end
            
            ScaledElement(~IdxNaN)  = Element;
            varargout(k)            = {ScaledElement};
        end
    end
end


