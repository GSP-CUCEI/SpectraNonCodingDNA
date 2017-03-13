% This script is to perform the step three of the normalization. The step is for each chromosome variable which already
% contain the Fourier coefficients of the RR for all the cell lines.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


load('C:\Users\Biomecanica1\Documents\MATLAB\WholeGenomeHumanRRSecond.mat', 'Chr1') % Variable with the Fourier coefficients of the genomic signal

for k=1:size(Chr1,1)
    for l=1:size(Chr1,2)      
        for m=1:size(Chr1(k,l).Matrix,2)
             
           Chr1AA(k,l).Matrix(m,:)=Chr1(k,l).Matrix(:,m)/max(Chr1(k,l).Matrix(:,m)); % Step three of the normalization
            
        end
        display(l)
    end
end

clear k l m Chr1
save HumanGenomeAA.mat Chr1AA '-v7.3'
clear Chr1AA


load('C:\Users\Biomecanica1\Documents\MATLAB\WholeGenomeHumanRRSecond.mat', 'Chr2') % Variable with the Fourier coefficients of the genomic signal

for k=1:size(Chr2,1)
    for l=1:size(Chr2,2)      
        for m=1:size(Chr2(k,l).Matrix,2)
             
           Chr2AA(k,l).Matrix(m,:)=Chr2(k,l).Matrix(:,m)/max(Chr2(k,l).Matrix(:,m)); % Step three of the normalization
             
        end
        display(l)
    end
end

clear k l m Chr2
save HumanGenomeAA.mat Chr2AA '-append' '-v7.3'
clear Chr2AA

load('C:\Users\Biomecanica1\Documents\MATLAB\WholeGenomeHumanRRSecond.mat', 'Chr3') % Variable with the Fourier coefficients of the genomic signal

for k=1:size(Chr3,1)
    for l=1:size(Chr3,2)      
        for m=1:size(Chr3(k,l).Matrix,2)
             
           Chr3AA(k,l).Matrix(m,:)=Chr3(k,l).Matrix(:,m)/max(Chr3(k,l).Matrix(:,m)); % Step three of the normalization
             
        end
        display(l)
    end
end

clear k l m Chr3
save HumanGenomeAA.mat Chr3AA '-append' '-v7.3'
clear Chr3AA

load('C:\Users\Biomecanica1\Documents\MATLAB\WholeGenomeHumanRRSecond.mat', 'Chr4') % Variable with the Fourier coefficients of the genomic signal

for k=1:size(Chr4,1)
    for l=1:size(Chr4,2)      
        for m=1:size(Chr4(k,l).Matrix,2)
             
           Chr4AA(k,l).Matrix(m,:)=Chr4(k,l).Matrix(:,m)/max(Chr4(k,l).Matrix(:,m)); % Step three of the normalization
            
        end
        display(l)
    end
end

clear k l m Chr4
save HumanGenomeAA.mat Chr4AA '-append' '-v7.3'
clear Chr4AA

load('C:\Users\Biomecanica1\Documents\MATLAB\WholeGenomeHumanRRSecond.mat', 'Chr5') % Variable with the Fourier coefficients of the genomic signal

for k=1:size(Chr5,1)
    for l=1:size(Chr5,2)      
        for m=1:size(Chr5(k,l).Matrix,2)
             
           Chr5AA(k,l).Matrix(m,:)=Chr5(k,l).Matrix(:,m)/max(Chr5(k,l).Matrix(:,m)); % Step three of the normalization
            
        end
        display(l)
    end
end

clear k l m Chr5
save HumanGenomeAA.mat Chr5AA '-append' '-v7.3'
clear Chr5AA

load('C:\Users\Biomecanica1\Documents\MATLAB\WholeGenomeHumanRRSecond.mat', 'Chr6') % Variable with the Fourier coefficients of the genomic signal

for k=1:size(Chr6,1)
    for l=1:size(Chr6,2)      
        for m=1:size(Chr6(k,l).Matrix,2)
             
           Chr6AA(k,l).Matrix(m,:)=Chr6(k,l).Matrix(:,m)/max(Chr6(k,l).Matrix(:,m)); % Step three of the normalization
            
        end
        display(l)
    end
end

clear k l m Chr6
save HumanGenomeAA.mat Chr6AA '-append' '-v7.3'
clear Chr6AA

load('C:\Users\Biomecanica1\Documents\MATLAB\WholeGenomeHumanRRSecond.mat', 'Chr7') % Variable with the Fourier coefficients of the genomic signal

for k=1:size(Chr7,1)
    for l=1:size(Chr7,2)      
        for m=1:size(Chr7(k,l).Matrix,2)
             
           Chr7AA(k,l).Matrix(m,:)=Chr7(k,l).Matrix(:,m)/max(Chr7(k,l).Matrix(:,m)); % Step three of the normalization
            
        end
        display(l)
    end
end

clear k l m Chr7
save HumanGenomeAA.mat Chr7AA '-append' '-v7.3'
clear Chr7AA

load('C:\Users\Biomecanica1\Documents\MATLAB\WholeGenomeHumanRRSecond.mat', 'Chr8') % Variable with the Fourier coefficients of the genomic signal

for k=1:size(Chr8,1)
    for l=1:size(Chr8,2)      
        for m=1:size(Chr8(k,l).Matrix,2)
             
           Chr8AA(k,l).Matrix(m,:)=Chr8(k,l).Matrix(:,m)/max(Chr8(k,l).Matrix(:,m)); % Step three of the normalization
            
        end
        display(l)
    end
end

clear k l m Chr8
save HumanGenomeAA.mat Chr8AA '-append' '-v7.3'
clear Chr8AA

load('C:\Users\Biomecanica1\Documents\MATLAB\WholeGenomeHumanRRSecond.mat', 'Chr9') % Variable with the Fourier coefficients of the genomic signal

for k=1:size(Chr9,1)
    for l=1:size(Chr9,2)      
        for m=1:size(Chr9(k,l).Matrix,2)
             
           Chr9AA(k,l).Matrix(m,:)=Chr9(k,l).Matrix(:,m)/max(Chr9(k,l).Matrix(:,m)); % Step three of the normalization
            
        end
        display(l)
    end
end

clear k l m Chr9
save HumanGenomeAA.mat Chr9AA '-append' '-v7.3'
clear Chr9AA

load('C:\Users\Biomecanica1\Documents\MATLAB\WholeGenomeHumanRRSecond.mat', 'Chr10') % Variable with the Fourier coefficients of the genomic signal

for k=1:size(Chr10,1)
    for l=1:size(Chr10,2)      
        for m=1:size(Chr10(k,l).Matrix,2)
             
           Chr10AA(k,l).Matrix(m,:)=Chr10(k,l).Matrix(:,m)/max(Chr10(k,l).Matrix(:,m)); % Step three of the normalization
            
        end
        display(l)
    end
end

clear k l m Chr10
save HumanGenomeAA.mat Chr10AA '-append' '-v7.3'
clear Chr10AA

load('C:\Users\Biomecanica1\Documents\MATLAB\WholeGenomeHumanRRSecond.mat', 'Chr11') % Variable with the Fourier coefficients of the genomic signal

for k=1:size(Chr11,1)
    for l=1:size(Chr11,2)      
        for m=1:size(Chr11(k,l).Matrix,2)
             
           Chr11AA(k,l).Matrix(m,:)=Chr11(k,l).Matrix(:,m)/max(Chr11(k,l).Matrix(:,m)); % Step three of the normalization
            
        end
        display(l)
    end
end

clear k l m Chr11
save HumanGenomeAA.mat Chr11AA '-append' '-v7.3'
clear Chr11AA

load('C:\Users\Biomecanica1\Documents\MATLAB\WholeGenomeHumanRRSecond.mat', 'Chr12') % Variable with the Fourier coefficients of the genomic signal

for k=1:size(Chr12,1)
    for l=1:size(Chr12,2)      
        for m=1:size(Chr12(k,l).Matrix,2)
             
           Chr12AA(k,l).Matrix(m,:)=Chr12(k,l).Matrix(:,m)/max(Chr12(k,l).Matrix(:,m)); % Step three of the normalization
            
        end
        display(l)
    end
end

clear k l m Chr12
save HumanGenomeAA.mat Chr12AA '-append' '-v7.3'
clear Chr12AA

load('C:\Users\Biomecanica1\Documents\MATLAB\WholeGenomeHumanRRSecond.mat', 'Chr13') % Variable with the Fourier coefficients of the genomic signal

for k=1:size(Chr13,1)
    for l=1:size(Chr13,2)      
        for m=1:size(Chr13(k,l).Matrix,2)
             
           Chr13AA(k,l).Matrix(m,:)=Chr13(k,l).Matrix(:,m)/max(Chr13(k,l).Matrix(:,m)); % Step three of the normalization
            
        end
        display(l)
    end
end

clear k l m Chr13
save HumanGenomeAA.mat Chr13AA '-append' '-v7.3'
clear Chr13AA

load('C:\Users\Biomecanica1\Documents\MATLAB\WholeGenomeHumanRRSecond.mat', 'Chr14') % Variable with the Fourier coefficients of the genomic signal

for k=1:size(Chr14,1)
    for l=1:size(Chr14,2)      
        for m=1:size(Chr14(k,l).Matrix,2)
             
           Chr14AA(k,l).Matrix(m,:)=Chr14(k,l).Matrix(:,m)/max(Chr14(k,l).Matrix(:,m)); % Step three of the normalization
            
        end
        display(l)
    end
end

clear k l m Chr14
save HumanGenomeAA.mat Chr14AA '-append' '-v7.3'
clear Chr14AA

load('C:\Users\Biomecanica1\Documents\MATLAB\WholeGenomeHumanRRSecond.mat', 'Chr15') % Variable with the Fourier coefficients of the genomic signal

for k=1:size(Chr15,1)
    for l=1:size(Chr15,2)      
        for m=1:size(Chr15(k,l).Matrix,2)
             
           Chr15AA(k,l).Matrix(m,:)=Chr15(k,l).Matrix(:,m)/max(Chr15(k,l).Matrix(:,m)); % Step three of the normalization
            
        end
        display(l)
    end
end

clear k l m Chr15
save HumanGenomeAA.mat Chr15AA '-append' '-v7.3'
clear Chr15AA

load('C:\Users\Biomecanica1\Documents\MATLAB\WholeGenomeHumanRRSecond.mat', 'Chr16') % Variable with the Fourier coefficients of the genomic signal

for k=1:size(Chr16,1)
    for l=1:size(Chr16,2)      
        for m=1:size(Chr16(k,l).Matrix,2)
             
           Chr16AA(k,l).Matrix(m,:)=Chr16(k,l).Matrix(:,m)/max(Chr16(k,l).Matrix(:,m)); % Step three of the normalization
            
        end
        display(l)
    end
end

clear k l m Chr16
save HumanGenomeAA.mat Chr16AA '-append' '-v7.3'
clear Chr16AA

load('C:\Users\Biomecanica1\Documents\MATLAB\WholeGenomeHumanRRSecond.mat', 'Chr17') % Variable with the Fourier coefficients of the genomic signal

for k=1:size(Chr17,1)
    for l=1:size(Chr17,2)      
        for m=1:size(Chr17(k,l).Matrix,2)
             
           Chr17AA(k,l).Matrix(m,:)=Chr17(k,l).Matrix(:,m)/max(Chr17(k,l).Matrix(:,m)); % Step three of the normalization
            
        end
        display(l)
    end
end

clear k l m Chr17
save HumanGenomeAA.mat Chr17AA '-append' '-v7.3'
clear Chr17AA

load('C:\Users\Biomecanica1\Documents\MATLAB\WholeGenomeHumanRRSecond.mat', 'Chr18') % Variable with the Fourier coefficients of the genomic signal

for k=1:size(Chr18,1)
    for l=1:size(Chr18,2)      
        for m=1:size(Chr18(k,l).Matrix,2)
             
           Chr18AA(k,l).Matrix(m,:)=Chr18(k,l).Matrix(:,m)/max(Chr18(k,l).Matrix(:,m)); % Step three of the normalization
            
        end
        display(l)
    end
end

clear k l m Chr18
save HumanGenomeAA.mat Chr18AA '-append' '-v7.3'
clear Chr18AA

load('C:\Users\Biomecanica1\Documents\MATLAB\WholeGenomeHumanRRSecond.mat', 'Chr19') % Variable with the Fourier coefficients of the genomic signal

for k=1:size(Chr19,1)
    for l=1:size(Chr19,2)      
        for m=1:size(Chr19(k,l).Matrix,2)
             
           Chr19AA(k,l).Matrix(m,:)=Chr19(k,l).Matrix(:,m)/max(Chr19(k,l).Matrix(:,m)); % Step three of the normalization
            
        end
        display(l)
    end
end

clear k l m Chr19
save HumanGenomeAA.mat Chr19AA '-append' '-v7.3'
clear Chr19AA

load('C:\Users\Biomecanica1\Documents\MATLAB\WholeGenomeHumanRRSecond.mat', 'Chr20') % Variable with the Fourier coefficients of the genomic signal

for k=1:size(Chr20,1)
    for l=1:size(Chr20,2)      
        for m=1:size(Chr20(k,l).Matrix,2)
             
           Chr20AA(k,l).Matrix(m,:)=Chr20(k,l).Matrix(:,m)/max(Chr20(k,l).Matrix(:,m)); % Step three of the normalization
            
        end
        display(l)
    end
end

clear k l m Chr20
save HumanGenomeAA.mat Chr20AA '-append' '-v7.3'
clear Chr20AA

load('C:\Users\Biomecanica1\Documents\MATLAB\WholeGenomeHumanRRSecond.mat', 'Chr21') % Variable with the Fourier coefficients of the genomic signal

for k=1:size(Chr21,1)
    for l=1:size(Chr21,2)      
        for m=1:size(Chr21(k,l).Matrix,2)
             
           Chr21AA(k,l).Matrix(m,:)=Chr21(k,l).Matrix(:,m)/max(Chr21(k,l).Matrix(:,m)); % Step three of the normalization
            
        end
        display(l)
    end
end

clear k l m Chr21
save HumanGenomeAA.mat Chr21AA '-append' '-v7.3'
clear Chr21AA

load('C:\Users\Biomecanica1\Documents\MATLAB\WholeGenomeHumanRRSecond.mat', 'Chr22') % Variable with the Fourier coefficients of the genomic signal

for k=1:size(Chr22,1)
    for l=1:size(Chr22,2)      
        for m=1:size(Chr22(k,l).Matrix,2)
             
           Chr22AA(k,l).Matrix(m,:)=Chr22(k,l).Matrix(:,m)/max(Chr22(k,l).Matrix(:,m)); % Step three of the normalization
            
        end
        display(l)
    end
end

clear k l m Chr22
save HumanGenomeAA.mat Chr22AA '-append' '-v7.3'
clear Chr22AA

load('C:\Users\Biomecanica1\Documents\MATLAB\WholeGenomeHumanRRSecond.mat', 'ChrX') % Variable with the Fourier coefficients of the genomic signal

for k=1:size(ChrX,1)
    for l=1:size(ChrX,2)      
        for m=1:size(ChrX(k,l).Matrix,2)
             
           ChrXAA(k,l).Matrix(m,:)=ChrX(k,l).Matrix(:,m)/max(ChrX(k,l).Matrix(:,m)); % Step three of the normalization
            
        end
        display(l)
    end
end

clear k l m ChrX
save HumanGenomeAA.mat ChrXAA '-append' '-v7.3'
clear ChrXAA
