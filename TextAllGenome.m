% In this script all the normalized coefficients matrixes for each RR, Cell Line, and Chromosome is converted to a 
% txt file.


HeadlinesRegulatorRegions={'ActPro','WeakPro','PoiPro','SEnhI','SEnhII','WeakEnhI','WeakEnhII',...
    'Ins','TXNTran','TXNElon','WeakTXN','Poly','Hetero','RepI','RepII'};
HeadlinesCellLines={'GM12878','H1esc','Hepg2','Hmec','Hsmm','Huevec','K562','Nhek','Nhlf'};

load('C:\Users\Biomecanica1\Documents\MATLAB\HumanGenomeAA.mat', 'Chr1AA')
 
for k=1:size(Chr1AA,1)
    for l=1:size(Chr1AA,2)
    
      A=Chr1AA(k,l).Matrix(2:end/2,:)';
      save(['Chr1AA' HeadlinesCellLines{k} HeadlinesRegulatorRegions{l} '.txt'],'A','-ascii')
      clear A
    end
end
clear Chr1AA

load('C:\Users\Biomecanica1\Documents\MATLAB\HumanGenomeAA.mat', 'Chr2AA')
 
for k=1:size(Chr2AA,1)
    for l=1:size(Chr2AA,2)
    
      A=Chr2AA(k,l).Matrix(2:end/2,:)';
      save(['Chr2AA' HeadlinesCellLines{k} HeadlinesRegulatorRegions{l} '.txt'],'A','-ascii')
      clear A
    end
end
clear Chr2AA

load('C:\Users\Biomecanica1\Documents\MATLAB\HumanGenomeAA.mat', 'Chr3AA')
 
for k=1:size(Chr3AA,1)
    for l=1:size(Chr3AA,2)
    
      A=Chr3AA(k,l).Matrix(2:end/2,:)';
      save(['Chr3AA' HeadlinesCellLines{k} HeadlinesRegulatorRegions{l} '.txt'],'A','-ascii')
      clear A
    end
end
clear Chr3AA

load('C:\Users\Biomecanica1\Documents\MATLAB\HumanGenomeAA.mat', 'Chr4AA')
 
for k=1:size(Chr4AA,1)
    for l=1:size(Chr4AA,2)
    
      A=Chr4AA(k,l).Matrix(2:end/2,:)';
      save(['Chr4AA' HeadlinesCellLines{k} HeadlinesRegulatorRegions{l} '.txt'],'A','-ascii')
      clear A
    end
end
clear Chr4AA

load('C:\Users\Biomecanica1\Documents\MATLAB\HumanGenomeAA.mat', 'Chr5AA')
 
for k=1:size(Chr5AA,1)
    for l=1:size(Chr5AA,2)
    
      A=Chr5AA(k,l).Matrix(2:end/2,:)';
      save(['Chr5AA' HeadlinesCellLines{k} HeadlinesRegulatorRegions{l} '.txt'],'A','-ascii')
      clear A
    end
end
clear Chr5AA

load('C:\Users\Biomecanica1\Documents\MATLAB\HumanGenomeAA.mat', 'Chr6AA')
 
for k=1:size(Chr6AA,1)
    for l=1:size(Chr6AA,2)
    
      A=Chr6AA(k,l).Matrix(2:end/2,:)';
      save(['Chr6AA' HeadlinesCellLines{k} HeadlinesRegulatorRegions{l} '.txt'],'A','-ascii')
      clear A
    end
end
clear Chr6AA

load('C:\Users\Biomecanica1\Documents\MATLAB\HumanGenomeAA.mat', 'Chr7AA')
 
for k=1:size(Chr7AA,1)
    for l=1:size(Chr7AA,2)
    
      A=Chr7AA(k,l).Matrix(2:end/2,:)';
      save(['Chr7AA' HeadlinesCellLines{k} HeadlinesRegulatorRegions{l} '.txt'],'A','-ascii')
      clear A
    end
end
clear Chr7AA

load('C:\Users\Biomecanica1\Documents\MATLAB\HumanGenomeAA.mat', 'Chr8AA')
 
for k=1:size(Chr8AA,1)
    for l=1:size(Chr8AA,2)
    
      A=Chr8AA(k,l).Matrix(2:end/2,:)';
      save(['Chr8AA' HeadlinesCellLines{k} HeadlinesRegulatorRegions{l} '.txt'],'A','-ascii')
      clear A
    end
end
clear Chr8AA

load('C:\Users\Biomecanica1\Documents\MATLAB\HumanGenomeAA.mat', 'Chr9AA')
 
for k=1:size(Chr9AA,1)
    for l=1:size(Chr9AA,2)
    
      A=Chr9AA(k,l).Matrix(2:end/2,:)';
      save(['Chr9AA' HeadlinesCellLines{k} HeadlinesRegulatorRegions{l} '.txt'],'A','-ascii')
      clear A
    end
end
clear Chr9AA

load('C:\Users\Biomecanica1\Documents\MATLAB\HumanGenomeAA.mat', 'Chr10AA')
 
for k=1:size(Chr10AA,1)
    for l=1:size(Chr10AA,2)
    
      A=Chr10AA(k,l).Matrix(2:end/2,:)';
      save(['Chr10AA' HeadlinesCellLines{k} HeadlinesRegulatorRegions{l} '.txt'],'A','-ascii')
      clear A
    end
end
clear Chr10AA

load('C:\Users\Biomecanica1\Documents\MATLAB\HumanGenomeAA.mat', 'Chr11AA')

for k=1:size(Chr11AA,1)
    for l=1:size(Chr11AA,2)
    
      A=Chr11AA(k,l).Matrix(2:end/2,:)';
      save(['Chr11AA' HeadlinesCellLines{k} HeadlinesRegulatorRegions{l} '.txt'],'A','-ascii')
      clear A
    end
end
clear Chr11AA

load('C:\Users\Biomecanica1\Documents\MATLAB\HumanGenomeAA.mat', 'Chr12AA')
 
for k=1:size(Chr12AA,1)
    for l=1:size(Chr12AA,2)
    
      A=Chr12AA(k,l).Matrix(2:end/2,:)';
      save(['Chr12AA' HeadlinesCellLines{k} HeadlinesRegulatorRegions{l} '.txt'],'A','-ascii')
      clear A
    end
end
clear Chr12AA

load('C:\Users\Biomecanica1\Documents\MATLAB\HumanGenomeAA.mat', 'Chr13AA')
 
for k=1:size(Chr13AA,1)
    for l=1:size(Chr13AA,2)
    
      A=Chr13AA(k,l).Matrix(2:end/2,:)';
      save(['Chr13AA' HeadlinesCellLines{k} HeadlinesRegulatorRegions{l} '.txt'],'A','-ascii')
      clear A
    end
end
clear Chr13AA

load('C:\Users\Biomecanica1\Documents\MATLAB\HumanGenomeAA.mat', 'Chr14AA')
 
for k=1:size(Chr14AA,1)
    for l=1:size(Chr14AA,2)
    
      A=Chr14AA(k,l).Matrix(2:end/2,:)';
      save(['Chr14AA' HeadlinesCellLines{k} HeadlinesRegulatorRegions{l} '.txt'],'A','-ascii')
      clear A
    end
end
clear Chr14AA

load('C:\Users\Biomecanica1\Documents\MATLAB\HumanGenomeAA.mat', 'Chr15AA')
 
for k=1:size(Chr15AA,1)
    for l=1:size(Chr15AA,2)
    
      A=Chr15AA(k,l).Matrix(2:end/2,:)';
      save(['Chr15AA' HeadlinesCellLines{k} HeadlinesRegulatorRegions{l} '.txt'],'A','-ascii')
      clear A
    end
end
clear Chr15AA

load('C:\Users\Biomecanica1\Documents\MATLAB\HumanGenomeAA.mat', 'Chr16AA')
 
for k=1:size(Chr16AA,1)
    for l=1:size(Chr16AA,2)
    
      A=Chr16AA(k,l).Matrix(2:end/2,:)';
      save(['Chr16AA' HeadlinesCellLines{k} HeadlinesRegulatorRegions{l} '.txt'],'A','-ascii')
      clear A
    end
end
clear Chr16AA

load('C:\Users\Biomecanica1\Documents\MATLAB\HumanGenomeAA.mat', 'Chr17AA')
 
for k=1:size(Chr17AA,1)
    for l=1:size(Chr17AA,2)
    
      A=Chr17AA(k,l).Matrix(2:end/2,:)';
      save(['Chr17AA' HeadlinesCellLines{k} HeadlinesRegulatorRegions{l} '.txt'],'A','-ascii')
      clear A
    end
end
clear Chr17AA

load('C:\Users\Biomecanica1\Documents\MATLAB\HumanGenomeAA.mat', 'Chr18AA')
 
for k=1:size(Chr18AA,1)
    for l=1:size(Chr18AA,2)
    
      A=Chr18AA(k,l).Matrix(2:end/2,:)';
      save(['Chr18AA' HeadlinesCellLines{k} HeadlinesRegulatorRegions{l} '.txt'],'A','-ascii')
      clear A
    end
end
clear Chr18AA

load('C:\Users\Biomecanica1\Documents\MATLAB\HumanGenomeAA.mat', 'Chr19AA')
 
for k=1:size(Chr19AA,1)
    for l=1:size(Chr19AA,2)
    
      A=Chr19AA(k,l).Matrix(2:end/2,:)';
      save(['Chr19AA' HeadlinesCellLines{k} HeadlinesRegulatorRegions{l} '.txt'],'A','-ascii')
      clear A
    end
end
clear Chr19AA

load('C:\Users\Biomecanica1\Documents\MATLAB\HumanGenomeAA.mat', 'Chr20AA')
 
for k=1:size(Chr20AA,1)
    for l=1:size(Chr20AA,2)
    
      A=Chr20AA(k,l).Matrix(2:end/2,:)';
      save(['Chr20AA' HeadlinesCellLines{k} HeadlinesRegulatorRegions{l} '.txt'],'A','-ascii')
      clear A
    end
end
clear Chr20AA

load('C:\Users\Biomecanica1\Documents\MATLAB\HumanGenomeAA.mat', 'Chr21AA')
 
for k=1:size(Chr21AA,1)
    for l=1:size(Chr21AA,2)
    
      A=Chr21AA(k,l).Matrix(2:end/2,:)';
      save(['Chr21AA' HeadlinesCellLines{k} HeadlinesRegulatorRegions{l} '.txt'],'A','-ascii')
      clear A
    end
end
clear Chr21AA

load('C:\Users\Biomecanica1\Documents\MATLAB\HumanGenomeAA.mat', 'Chr22AA')
 
for k=1:size(Chr22AA,1)
    for l=1:size(Chr22AA,2)
    
      A=Chr22AA(k,l).Matrix(2:end/2,:)';
      save(['Chr22AA' HeadlinesCellLines{k} HeadlinesRegulatorRegions{l} '.txt'],'A','-ascii')
      clear A
    end
end
clear Chr22AA

load('C:\Users\Biomecanica1\Documents\MATLAB\HumanGenomeAA.mat', 'ChrXAA')
 
for k=1:size(ChrXAA,1)
    for l=1:size(ChrXAA,2)
    
      A=ChrXAA(k,l).Matrix(2:end/2,:)';
      save(['ChrXAA' HeadlinesCellLines{k} HeadlinesRegulatorRegions{l} '.txt'],'A','-ascii')
      clear A
    end
end
clear ChrXAA
