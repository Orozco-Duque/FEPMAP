function [ Th1_V,NumP ] = AdaptativeThreshold2(signal,N,Th1)

   
    NPKI=0.2*Th1;
    SPKI=1.2*Th1;
    NumP=0;
    for n1=1:length(signal)-N
        PEAKI=max(signal(n1:n1+N));
        if PEAKI>Th1
            SPKI=0.125*PEAKI+0.875*SPKI;
            %Th1=NPKI+0.25*(SPKI-NPKI);
            Th1=NPKI+0.35*(SPKI-NPKI);
            NumP=NumP+1;
        else
        NPKI=0.125*PEAKI+0.875*NPKI;
        end
        %Th1=Th1-0.01*PEAKI;
        Th1_V(n1)=Th1;
        Th1_V(length(Th1_V):length(signal))=Th1_V(length(Th1_V));
    end
 end

