
%% ofdm_tx_longpn
function ofdm_tx_longpn(pn_len, config)
    % clear all;

    if nargin < 1, pn_len = 128; end
    if nargin < 2, config = 1; end


    output_dir = './tx_sound/';

    % paramters
    % config=1;
    if config==1
        Nfft=128;
        Ncp=32;
        symbolRate=1/2205;
        Ts=1/44100;
        Fs=1/Ts;
        fc=18000;
        rollOff=0.25;
        nSamp=20;
        filterSpan=10;
        % preambleSeq=m_sequence([1 0 0 0 0 1 1 1]);
        preambleSeq=m_sequence([1 0 0 0 0 1 1 1 0 1 0 1 0 0 0 1 0]);
    elseif config==2
        Nfft=128;
        Ncp=32;
        symbolRate=1/2205;
        Ts=1/44100;
        Fs=1/Ts;
        fc=16000;
        rollOff=0.25;
        nSamp=20;
        filterSpan=10;
        % preambleSeq=m_sequence([1 0 0 0 0 1 1 1]);
        preambleSeq=m_sequence([1 0 0 0 0 1 1 1 0 1 0 1 0 0 0 1 0]);
    else
        error('wrong config');
    end



    % generate data
    data=preambleSeq;
    % data=data(1:Nfft);
    data = data(1:min(end,pn_len));
    data = reshape(data(1:Nfft*floor(length(data)/Nfft)), Nfft, []);
    % modulate with BPSK
    data=data*2-1;

    % generate OFDM symbol
    data_fft = ifft(data,Nfft);
    data_fft = reshape(data_fft, 1, []);
    preamble_size = length(data_fft);
    fprintf(' Size of Preamble: %d (~%fs)\n', preamble_size, preamble_size*Nfft*Ts);
    % plot(abs(data_fft));
    ofdmSymbol=[data_fft(Nfft-Ncp+1:Nfft), data_fft];
    plot(abs(ofdmSymbol));

    % pulse shaping
    firConfig=fdesign.pulseshaping(nSamp,'Raised Cosine','Nsym,Beta',filterSpan,rollOff);
    fir=design(firConfig);
    firValue=fir.Numerator;
    preamble=upfirdn(ofdmSymbol, firValue, nSamp);
    analogData=[zeros(1,Fs/100) preamble zeros(1,Fs/100)];
    plot(abs(preamble));
    % disp(numel(analogData));
    % return

    fprintf('  one packet size = %dx%d\n', size(analogData));

    % generating multiple symbols
    % intervalSignal=zeros(1,ceil(Fs/3));
    % intervalSignal=zeros(1, length(analogData));
    intervalSignal=zeros(1, ceil(Fs*2));
    fprintf('  one packet size with interval = %dx%d\n', size([intervalSignal, analogData, intervalSignal]));
    Ns=60;
    oneSymbol = [intervalSignal, analogData, intervalSignal]; 
    multiSymbol=[];
    for i=1:Ns
        multiSymbol=[multiSymbol, intervalSignal, analogData, intervalSignal]; 
        % fprintf('  one packet size with interval = %dx%d\n', size([intervalSignal, analogData, intervalSignal]));
    end
    analogData=multiSymbol;

    % upconverting
    T=numel(analogData);
    analogData=real(analogData.*exp(1i*2*pi*fc*(1:T)*Ts));
    analog_one_symbol = real(oneSymbol.*exp(1i*2*pi*fc*(1:numel(oneSymbol))*Ts));
    % dataSpectrum=fft(analogData);
    % dataSpectrum=circshift(dataSpectrum,[0 floor(T/2)]);
    % Th=(T-1)/2;
    % plot((-Th:Th)/Th*Fs/2,abs(dataSpectrum));

    % scaling
    analogData=analogData/max(abs(analogData))*0.9;
    analog_one_symbol = analog_one_symbol / max(abs(analog_one_symbol)) * 0.9;
    plot(abs(analogData));

    % % model channel effects
    % htap=5;
    % h=zeros(1,htap*nSamp);
    % h(1:nSamp:(htap-1)*nSamp+1)=[1 5 2 1 1];
    % analogData=conv(analogData,h);
    % plot(abs(analogData));

    % write the results
    save([output_dir 'preamble.' num2str(config) '.' num2str(fc) '.' num2str(preamble_size) '.mat'],'preamble');
    save([output_dir 'txch_signal.' num2str(config) '.' num2str(fc) '.' num2str(preamble_size) '.mat'],'analogData');
    save([output_dir 'txch_1signal.' num2str(config) '.' num2str(fc) '.' num2str(preamble_size) '.mat'],'analog_one_symbol');


    % generate audio file
    signal_d=zeros(T,2);
    signal_d(:,2)=analogData.';
    audiowrite([output_dir 'ofdm.' num2str(config) '.' num2str(fc) '.' num2str(preamble_size) '.wav'], signal_d, Fs, 'BitsPerSample', 16);
end