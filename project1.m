waveFile='test.wav';

% read file
[y, fs] = audioread(waveFile);

% initail constant
frameTime = 20; % (ms)
overlapTime = 10; % (ms)
silenceTime = 100; % (ms)
wantedAutocorrIndex = 400;
maxFreq = 200;
minFreq = 50;

frameSize = round(fs * frameTime / 1000);
shiftSize = round(fs * (frameTime - overlapTime) / 1000);
silenceSize = round(fs * silenceTime / 1000);

frameNum = ceil((length(y) - frameSize) / shiftSize) + 1;
silenceNum = ceil((silenceSize - frameSize) / shiftSize + 1);

% initial window function
win = window(@hamming, frameSize);

% initial sgn to compute ZCR
sgn = sign(y);
sgn = sgn + (sgn == 0);

% compute Energy, ZCR, pitch
energy = zeros(1, frameNum);
zcr = zeros(1, frameNum);
pitch = zeros(1, frameNum);

for i = 1 : frameNum
    startIndex = (i - 1) * shiftSize + 1;
    endIndex = min(startIndex + frameSize - 1, length(y) -  1);
    
    % hamming window
    tmpEnergy = zeros(1, frameSize);
    tmpZcr = abs(sgn(startIndex + 1 : endIndex + 1) - sgn(startIndex : endIndex));
    for j = startIndex : endIndex
        idx = j - startIndex + 1;
        tmpEnergy(idx) = y(j) * win(idx);
        tmpZcr(idx) = tmpZcr(idx) * win(idx);
    end
    
    energy(i) = sum(tmpEnergy .^ 2);
    zcr(i) = sum(tmpZcr);
    
    % compute autocorrelation
    autocorrelation = zeros(1, frameSize);
    for k = 0 : frameSize
        tmpAutocorr = zeros(1, frameSize);
        for j = 0 : frameSize - k
            if startIndex + j + k > endIndex
                break;
            end
        tmpAutocorr(j + 1) = y(startIndex + j) * y(startIndex + j + k);
        end
        autocorrelation(k + 1) = sum(tmpAutocorr);
    end
    if i == wantedAutocorrIndex
        wantedAutocorr = autocorrelation;
    end
    
    [pks, locs] = findpeaks(autocorrelation(1 + round(fs / maxFreq) : end));
    if isempty(locs) || energy(i) < 1
        pitch(i) = 0;
    else
        pitch(i) = fs / (locs(find(pks == max(pks), 1)) + round(fs / maxFreq));
    end
end

% compute ITU, ITL
IMX = max(energy);
IMN = mean(energy(1:silenceNum));
I1 = 0.003 * (IMX - IMN) + IMN;
I2 = 4 * IMN;
ITL = min(I1, I2);
ITU = 5 * ITL;

% compute IZCT
IF = (25 / 10 * frameTime) * 2; % 25crossing / 10ms
IZC = mean(zcr(1:silenceNum));
stdIZC = std(zcr(1:silenceNum));
IZCT = min(IF, IZC + 2 * stdIZC);

% end-point detection
N1 = find(energy >= ITU, 1);
tmp = find(energy(1 : N1) <= ITL, 1, 'last');
if ~isempty(tmp)
    N1 = tmp;
end
tmp = find(zcr(1 : N1) >= IZCT, 1, 'last');
if ~isempty(tmp)
    N1 = tmp;
end

N2 = find(energy >= ITU, 1, 'last');
tmp = find(energy(N2 : end) <= ITL, 1);
if ~isempty(tmp)
    N2 = N2 - 1 + tmp;
end
tmp = find(zcr(N2 : end) >= IZCT, 1);
if ~isempty(tmp)
    N2 = N2 - 1 + tmp;
end

% plot
frameAxis = 1 : frameNum;
timeAxis = 1 : (frameNum - 1) / (length(y) - 1) : frameNum;

subplot(6, 1, 1);
plot(timeAxis, y);
title('waveform');

subplot(6, 1, 2);
plot(frameAxis, energy);
title('energy');

subplot(6, 1, 3);
plot(frameAxis, zcr);
title('ZCR');

subplot(6, 1, 4);
plot(wantedAutocorr);
titleStr = sprintf('autocorrelation of frame%d', wantedAutocorrIndex);
title(titleStr);

subplot(6, 1, 5);
plot(frameAxis, pitch);
title('pitch');

subplot(6, 1, 6);
plot(timeAxis, y);
line([N1, N1], ylim, 'Color', 'red');
text(N1,  max(ylim), 'Beginning', 'Color', 'red', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
line([N2, N2], ylim, 'Color', 'green');
text(N2, max(ylim), 'Ending', 'Color', 'green', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
