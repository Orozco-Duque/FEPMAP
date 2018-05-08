%%CFEmean2
%Compute the proportion of each EGM signal holding fractionated electrical 
%activity, and is calculated by detection of deflection times
%Written by: Andrés Orozco - UPB - Medellín, Colombia. 2013
%Input:
%-signal
%-sensitivity: threshold for the mimimum peak to peak amplitud in a deflection 
%-windows_width: maximum value for delflection intervale width (to avoid
%ventricular potentials)
%-refractory_period: Minimum value between deflections 
%-fig: 1 to plot the result - 0 not plotting.
%Output
%-IF: Auxiliar signal with non-zero values in segments with fragmented
%activity
%-CFEm: mean of fragmented intervales duration

function [y CFEm width anchoy Act IF]=FragmentedAct(OriginalSignals,sensitivity,window_width,refractory_period,fig)

if nargin < 5

	fig=0;
end
if nargin < 4

	refractory_period=30;
end
if nargin < 3
	window_width=20;

end
if nargin < 2
	sensitivity=0.1*max(OriginalSignals);
end


if size(OriginalSignals,1)>size(OriginalSignals,2)
    OriginalSignals=OriginalSignals';
end

    signal=OriginalSignals;
   % signal1=signal1/max(abs(signal1));
    
     
    N=1:length(signal);
    difsignal=-diff(signal);
    iniD=[];
    endD=[];
    iniDm=[];
    endDm=[];
    flagD=0;
    flagDm=0;
     
    y=zeros(1,length(signal));
    
    %Detection of maximum-minimum intervale
    for i=1:length(signal)-1
        if difsignal(i)>0 && flagD==0
            iniD=[iniD i];
            flagD=1;
        end
        if difsignal(i)<0 && flagD==1
            endD=[endD i];
            flagD=0;
        end
    end
    
    %Detecting mimimum-maximum intervale
    
    for i=1:length(signal)-1
        if difsignal(i)>0 && flagDm==1
            iniDm=[iniDm i];
            flagDm=0;
        end
        if difsignal(i)<0 && flagDm==0
            endDm=[endDm i];
            flagDm=1;
        end
    end
    
    
    
    if length(iniD)~=length(endD) % Maximum-minimum intervals 
        %and minimum-maximum must have the same length
        len=min(length(iniD),length(endD));
        iniD=iniD(1:len);
        endD=endD(1:len);
    end
    
    
    
    if length(iniDm)~=length(endDm) 
        len=min(length(iniDm),length(endDm));
        iniDm=iniDm(1:len);
        endDm=endDm(1:len);
    end
   
    %For each maximum-minimum intervale
    j=1;
    for i=1:length(iniD)
        if (endD(i)-iniD(i))<window_width
            endD2(j)=endD(i);
            iniD2(j)=iniD(i);
            j=j+1;
        end
    end
    
    %For each mimimum-maximum intervale
    j=1;
    for i=1:length(iniDm)
        if (iniDm(i)-endDm(i))<window_width
            endD2m(j)=endDm(i);
            iniD2m(j)=iniDm(i);
            j=j+1;
        end
    end
    

    j=1;
    s2=abs(signal(iniD2)-signal(endD2));
    for i=1:length(s2)
        if s2(i)>sensitivity
            iniD3(j)=iniD2(i);
            endD3(j)=endD2(i);
            s3(j)=s2(i);
            j=j+1;
        end
    end
    
    

    j=1;
    try
    s2m=abs(signal(endD2m)-signal(iniD2m));
    
        
    for i=1:length(s2m)
        if s2m(i)>sensitivity
            iniD3m(j)=iniD2m(i);
            endD3m(j)=endD2m(i);
            s3m(j)=s2m(i);
            j=j+1;
        end
    end
    catch
    iniD3m=0;
    endD3m=0;
    s2m=0;
    end
    try
    iniD4(1)=iniD3(1);
     picos=length(iniD3m);
    picos=length(iniD3);
    
    catch
        CFEm=0;
        IF=0;
        AmpMean=0;
        Act=[];
        picos=0;
        width=0;
        anchoy=0;
        return;
    end
    
    
    
    
    iniD=iniD3;
    endD=endD3;
    
   
    for i=1:length(endD)
        
        
        y(iniD(i):endD(i))=s3(i);  %auxiliar signals with values different 
        %to zero in segments with maximum-minimum intervals        
    end
    
    if exist('s3m')
        for i=1:length(endD3m)
       
            y(endD3m(i):iniD3m(i))=s3m(i);  
        
        end
    end
    
    %Code to put together segments 
    
    flagy=0;
    timeIseg=0;
    for i=1:length(y)
        if y(i)>0 && flagy==0
            flagy=1;
        end
        if y(i)==0 && flagy==1
            
            flagy=2;
            iniseg=i;
            yant=y(i-1);
            timeIseg=0;
        end
        if flagy==2
            timeIseg=timeIseg+1;
        end
        if flagy==2 && y(i)>0
            flagy=0;
            if timeIseg<refractory_period
                y(iniseg:i)=min([yant y(i)]);
                
            end
        end
    end
            
    
    
    %Detection of values less than 20% of the maximum peak
    
    
    flagy=0;
    for i=1:length(y)
        if y(i)>0 && flagy==0
            iniseg=i;
            flagy=1;
        end
        if y(i)==0 && flagy==1
            endseg=i;
            flagy=0;
            yseg=y(iniseg:endseg);
            maxSeg=max(yseg);
            yseg(yseg<(maxSeg*0.15))=0;
            y(iniseg:endseg)=yseg;
        end
    end
    
    %Put together
    
    flagy=0;
    timeIseg=0;
    for i=1:length(y)
        if y(i)>0 && flagy==0
            flagy=1;
        end
        if y(i)==0 && flagy==1
            
            flagy=2;
            iniseg=i;
            yant=y(i-1);
            timeIseg=0;
        end
        if flagy==2
            timeIseg=timeIseg+1;
        end
        if flagy==2 && y(i)>0
            flagy=0;
            if timeIseg<refractory_period
                y(iniseg:i)=min([yant y(i)]);
                
            end
        end
    end
    
    flagy=0;
   
    j=0;
     for i=1:length(y)
         if y(i)>0 && flagy==0
             j=j+1;
             flagy=1;
             anchoy(j)=0;
         end
         if flagy==1
             anchoy(j)=anchoy(j)+1;
         end
         if y(i)==0 && flagy==1
             flagy=0;
         end
     end
         
           
    
    width=length(y(y>0));
    
    
    %Validate the maxima
    j=1;
    for i=1:length(iniD)
        Ac=round(iniD(i)+(endD(i)-iniD(i))/2);
        if y(Ac)>0
            Act(j)=Ac;
            j=j+1;
        end
    end
        
    if fig==1
        plot(N,signal,N,y,'k--',N(Act),signal(Act),'r+')
        
    end
    
    for i=1:length(Act)-1
        IF(i)=Act(i+1)-Act(i);
        
    end
    
    
    try
    CFEm=mean(IF);
    catch
        CFEm=0;
        IF=0;
    end
end
