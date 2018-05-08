
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%Function EGMFeatures4 compute morphological and non-linear features from a
%single EMG signal.
%
%Inputs:
%Signal: A vector that contain the signal
%Fs: sample frequency in Hz
%FL: Low cutoff frequency
%FH: High cutoff frequency
%th: Threshold to be used in the activity detection
%
%Outputs:
%FAr: Fractionated activity
%SI: Similarity index
%MF: Multifractal h-fluctuation index
%ApEn: Approximate entropy
%
%Written by Andres Orozco-Duque
%Biomedical Engineering Lab
%Instituto Tecnológico Metropolitano - Medellín - Colombia

function [FAr SI MF ApEn ]=EGMfeatures4(signal,fs,th)
    
if nargin < 3
    %If threshold is not available, the default is 0.1V.
	th=0.1;
end

    %%
    %Computation of filter parameters for Local Activation Waves (LAW) detection
    fcbp=[30 300];
    fclp=20;
    orden=4;
    wn1=fcbp/(fs/2);                        %Frecuencias normalizadas pasabanda
    wn2=fclp/(fs/2);                        %Frecuencia normalizada pasabajo
    d1=fdesign.bandpass('N,F3dB1,F3dB2',orden,wn1(1,1),wn1(1,2));   % Butter
    d4=fdesign.lowpass('N,F3db',orden,wn2);     % Butter, maxflat
    BPcoef=design(d1,'butter');
    LPcoef=design(d4,'butter'); 
    
        %%
        %Fractionated activity describes the proportion of each EGM signal 
        %presenting continuous electrical activity.
        [y]=CFEmean2(signal,th,20,50,0);
        FAr=length(find(y~=0))/length(y);   %radio de actividad fragmentada
                   
        %%        
        %Similarity index
         [Act_loc,sensitivityT,eAtime]=dectectLAW(signal,fs,0,BPcoef,LPcoef,th);
        if ~isempty(Act_loc)
           [SI] =similarity_fun2(signal,Act_loc);
        else
           SI=0;
        end
        
        %%
        %approximate entropy
        ApEn=appen(signal,3,0.38*std(signal));
       
        %%
        %Multifractal h-fluctuation index
        m=2;
        scres=40;
        scmax=256;
        scmin=64;
        exponents=linspace(log2(scmin),log2(scmax),scres);
        scale=round(2.^exponents);
        q=[-3:0.1:3];
        [Hq,tq,hq,Dq,Fq]=MFDFA1(signal,scale,q,m,0);
         [~,MF,~,~,~,~,~]=...
                         multifractalFeatures(hq,Dq,30,0);
         

end

