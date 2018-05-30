%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%DetectLAW: Function for LAW detection 
%
% Input
% signalP:
% fs:
% plotting

function [Act_loc1,th,eAtime] = dectectLAW(signalP,fs,plotting,BPcoef,LPcoef,threshold)


if nargin < 3

	plotting=1;
end

if nargin < 2

	fs=2034.5;
end

if size(signalP,2)==1
    signalP=signalP';
end

f1=filter(BPcoef,signalP); %bandpass filter
rec=abs(f1);                %rectification
f2=filter(LPcoef,rec);    %lowpass filter to construct the envelope


[th,NumP ] = AdaptativeThreshold2(f2,140,threshold/2);
        
        
        


if size(th,2)==1
    th=th';
end
wind=0.050*fs;
wind=round(wind);
Act_loc=[];
poss_max=(f2>th); % It detects each point over the threshold
leftmax = find(diff([0 poss_max])==1);% to the left
rigthmax= find(diff([poss_max 0])==-1);% to the righ
for i=1:length(leftmax)
    if leftmax(i)<=wind || rigthmax(i)>=length(signalP)-wind
         aux=signalP(leftmax(i):rigthmax(i));
         posaux=find(aux==max(signalP(leftmax(i):rigthmax(i))));
         Act_loc(i) = posaux(1);
         Act_loc(i) = Act_loc(i)-1+leftmax(i); % add offset 
    else
        aux=signalP(leftmax(i)-wind:rigthmax(i)+wind);
        posaux=find(aux==max(signalP(leftmax(i)-wind:rigthmax(i)+wind)));
        Act_loc(i) = posaux(1);
        Act_loc(i) = Act_loc(i)-1+leftmax(i)-wind; % add offset 
    end
end

eAtime=sum(rigthmax-leftmax)./fs;


a=0;
Act_loc1=Act_loc; 



while ~isempty(find(diff(Act_loc1)<round(0.07*fs))) % It delect minimal points with separation less than 70 ms 
for i=1:length(Act_loc1)-1
    if (Act_loc1(i+1)-Act_loc1(i))<round(0.07*fs)
        aux=([signalP(Act_loc1(i)) signalP(Act_loc1(i+1))]);
        loc=find(aux==min([signalP(Act_loc1(i)) signalP(Act_loc1(i+1))]));
        a=a+1;
        borr(a)=i+loc(1)-1;
    end
         
end 
Act_loc1(borr)=[]; % Activation points
a=0;
borr=[];

end



if plotting==1
    t=[1/fs:1/fs:length(signalP(1,:))/fs];
    subplot(2,1,1),plot(t,f2),hold on,plot(t,th,'r-.')
    subplot(2,1,2),plot(t,signalP,t(Act_loc1),signalP(Act_loc1),'ro')
 % plot(t,signalP,t(Act_loc1),signalP(Act_loc1),'ro',t,f2,'--',t,th,'r-.')
end
end