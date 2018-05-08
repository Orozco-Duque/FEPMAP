%multifractalFeatures
%
%Function to generate a feature with information about the behavior of h(q)
%in multifractal analisys
%
%Inputs:
%hq: h(q) is the hq result of multifractal analisys by DFA or WTMM
%Dq: D(q) 
%lwindow: the length of the analisys windows related to q
%standar: if 1, spectrum is standarize
%
%Outputs:
%d2hq: second derivative of h(q) 
%mfeature: sum of the desviation of d2hq to zero (multifractal feature)
%hqwidth: h(q=-inf)-h(q=+inf) note: inf is the last of firts position
%poshqmax: Desviation in the position of Dq(max) respect to q=0;
%hqmax: hq(q=Dq(max))
%W80: h2-h1 / D(h1)=0.8D(hmax) and D(h2)=0.8D(hmax)
%dhleft: h(q=Dq(max))-hq(q=inf)
%
%
%Implemented by:
%Andrés Orozco-Duque
%Universidad Pontificia Bolivariana - Centro de Bioingeniería
%Medellin-Colombia 2013
%
%Supported by:
%Daniel Novak
%Czech Technical University in Prague

function [d2hq mfeature hqwidth poshqmax hqmax W80 dhleft dharea]=multifractalFeatures(hq,Dq,lwindow,standar)

    if nargin < 4
       standar=1;
    end
    if nargin < 3
       lwindow=length(hq);
    end

    if standar==1
        hq=(hq-mean(hq))/std(hq);
    end
    lhalf=floor(lwindow/2);
 %   h1n(:,i)=Dn{i}(:,1);
  %  hq0delta(i)=h1(100,i)-h1n(100,i);;
    pm=floor(length(hq)/2);  %h midpoint
    poshqmax=find(Dq==max(Dq));
    poshqmax=poshqmax(1);
    hqmax=hq(poshqmax);
    h2=hq(max(find(Dq(1:poshqmax)<(0.8*Dq(poshqmax)))));
    if isempty(h2)
        h2=hq(1);
    end
    h1=hq(min(find(Dq(poshqmax:end)<(0.8*Dq(poshqmax))))+poshqmax-1);
    if isempty(h1)
        h1=hq(end);
    end
    W80=h2-h1;
    hq0=hq(pm);
    hqright=hq(1);
    hqleft=hq(end);
    hqwidth=hq(1)-hq(end);
    dhright=(hqright-hqmax)/hqwidth;
    dhleft=(hqmax-hqleft)/hqwidth;
    dharea=sum(Dq);
 
    
    
    hwindow=hq(pm-lhalf+1:pm+lhalf); %I make a window
    d2hq=diff(diff(hwindow));
    mfeature=sum((d2hq).^2)/length(d2hq);
    
    poshqmax=poshqmax-pm;
  
end
