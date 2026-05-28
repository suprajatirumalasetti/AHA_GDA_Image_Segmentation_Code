function [Dt] = M_CEM(I_energy,Th)
    Th = sort(Th);
      I_energy = double(I_energy);  
hn = I_energy ./ (sum(I_energy) + eps);
    for ii=1:length(Th)+1
        m1=0;
        m2=0;
        if ii==1
            v1=1;
            v2=Th(ii)-1;
        elseif ii==(length(Th)+1)
            v1=Th(ii-1);
            v2=256;
        else
            v1=Th(ii-1);
            v2=Th(ii)-1;
        end
        v1 = max(1, v1);
        v2 = min(256, v2);
        
        for i = v1:v2
            if i >= 1 && i <= length(hn)
                m1 = m1 + (round(i) * hn(round(i)));
                m2 = m2 + hn(round(i));

            end
        end
        if m2 == 0
            miu(ii)=m1/(m2+eps);
        else
            miu(ii)=m1/m2;
        end
        Entro = double(0);
        for i = v1:v2 
            if i >= 1 && i <= 256
                if miu(ii) == 0
                   Entro = double(Entro + round(i) * hn(round(i)) * log(miu(ii) + eps));
                else
                     Entro = double(Entro + (round(i) * hn(round(i)) * log(miu(ii))));
                end
            end
        end
        Entropy(ii) = Entro;     
    end
        imEntropy = double(0);
    for i = 1:256
        if i >= 1 && i <= 256
            imEntropy = double(imEntropy+(i*hn(i)*log(i)));
        end
    end
    Entemp=0;
    for i =1: length(Entropy)
        Entemp = Entemp + Entropy(i);
    end
    Dt = imEntropy - Entemp;
end