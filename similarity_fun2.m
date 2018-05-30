% A Method for Quantifying Atrial
% Fibrillation Organization Based
% on Wave-Morphology Similarity

function [ro, s1 pa pb x] =similarity_fun2(signal1,Act_loc1)
f=0;
epsilon=0.8;
 f=f+1;
wind=90;% To generate a window 
wind2=floor(wind/2); % the middle point of the window

    signal=signal1;
    
   % [Act_loc1] = dectect_activsim(signal,0,3,1,'mexh',150,20); % entrega las activaciones

    if length(Act_loc1)>5
    Act_loc1(end)=[]; % se elimina la ultima muestra de activacion por si no está completa y no dañe la comparacion
    Act_loc1(1)=[]; % se elimina la primera muestra de activacion por si no está completa y no dañe la comparacion
    
    elseif length(Act_loc1)>3 && length(Act_loc1)<=5
     Act_loc1(1)=[];
    end
    
    N=length(Act_loc1);
    x=zeros(N,wind);
    s=zeros(N,wind);
  
    if Act_loc1(1)>wind2 && Act_loc1(1)<(length(signal)-wind2)% 
    x(1,:)=signal(-(wind2-1)+Act_loc1(1):Act_loc1(1)+wind2); % 
    s(1,:)=x(1,:)/sqrt(sum(x(1,:).^2)); % n
    else
    x(1,:)=signal(1:wind); % 
    s(1,:)=x(1,:)/sqrt(sum(x(1,:).^2));
    end
    
    if Act_loc1(end)<(length(signal)-wind2) % 
    x(N,:)=signal(-(wind2-1)+Act_loc1(N):Act_loc1(N)+wind2);
    s(N,:)=x(N,:)/sqrt(sum(x(N,:).^2));
    else
    x(N,:)=signal(end-wind+1:end);
    s(N,:)=x(N,:)/sqrt(sum(x(N,:).^2));
    end
    
    for i=2:N-1
    x(i,:)=signal(-(wind2-1)+Act_loc1(i):Act_loc1(i)+wind2); % extrae los segmentos sin incluir el primero ni el ultimo
    s(i,:)=x(i,:)/sqrt(sum(x(i,:).^2));
    end
     
    
    a=0;
    pa=0;
    pb=0;
    p1=0;
    s1=[];
    for i=1:size(s,1)
        for j=i+1:size(s,1)
           a=a+1;
           [c,lags]=xcorr(s(i,:),s(j,:));
           [~, loc]=max(c);
           desf=lags(loc);
           x1(i,:)=x(i,:);
           s1(i,:)=x1(i,:)/sqrt(sum(x1(i,:).^2));
           
           if Act_loc1(j)+wind2-desf<length(signal)&& Act_loc1(j)-wind2-desf>1
           x1(j,:)=signal(-(wind2-1)+Act_loc1(j)-desf:Act_loc1(j)+wind2-desf); 
           s1(j,:)=x1(j,:)/sqrt(sum(x1(j,:).^2));
          
           p1(a)=1-real(acos(s1(i,:)*s1(j,:)'));  % 
           pa(a)=heaviside(epsilon-real(acos(s1(i,:)*s1(j,:)'))); 
           
           corre=corrcoef(s1(i,:),s1(j,:)); % 
           pb(a)=abs(corre(1,2));
           
           else
           pa(a)=0;
           pb(a)=0;
            end
            
        end
    end
  
 ro=[mean(pa)];% mean(pb)];% mean(pc) mean(pd)];

end



