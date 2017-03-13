% To create the variable used in this script we downloaded the archives from ENCODE database related 
% to ChromHMMseg for all the cell lines. Also downloaded the fasta archives used as reference in the performance
% of the ChromHMMseg algorithm reported in Ernst et al. (2011). 
% We extracted the sequences of previously mentioned build of the human genome with the indexes of
% the ChromHMMseg files.

% Paths in this script are the used in this work.

% The next code was performed for each chromosome (1 to 22 and X chromosome). 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This section is to call all the variable needed
% Cell line B-Lymphoblastoid cells

counter=zeros(23,15,9);
load('C:\Users\Biomecanica1\Documents\MATLAB\Final\SizeRR.mat', 'SizeGM12878') % Variable with the length of the RR genomic signals
load('C:\Users\Biomecanica1\Documents\MATLAB\DNA\GM12878.mat', 'GM12878') % Variable with the extracted RR DNA sequences
load('C:\Users\Biomecanica1\Documents\MATLAB\Final\WholeGenomeStatistic.mat', 'WholeMedian') % This is the median of the whole genomic signals
save WholeGenomeHumanRR.mat counter '-v7.3'
save WholeGenomeHumanRRE.mat counter 

% In this section, the DNA sequence are mapped to a genomic signal and also the step one and two of the normalization
% is performed.

for l=1:15
   for m=1:size(GM12878(1).Region(l).Zone,2)

       R=find(upper(GM12878(1).Region(l).Zone(m).Sequence)=='N'); % To evaluate if a sequence has an any nucleotide within
       counter(1,l,1)=0;
           
       if (isempty(R)) 
               
           Signal=mapping(GM12878(1).Region(l).Zone(m).Sequence,3); % Neighbor joining mapping
           Fourier=abs(fft(Signal)); % FFT calculation
           Chr1(1,l).Matrix(m-counter(1,l,1),:)=interp1(1:size(Fourier,2),Fourier,linspace(1,size(Fourier,2),WholeMedian)); % Step one of the normalization
           Chr1(1,l).Matrix(m-counter(1,l,1),:)=Chr1(1,l).Matrix(m-counter(1,l,1),:)/SizeGM12878{1,l}(m); % Step two of the normalization
               
       else
               
           counter(1,l,1)=counter(1,l,1)+1; % To count if a sequence has an any nucleotide within
           
       end
           
       clear Signal Fourier R

   end
   display(l)
   save WholeGenomeHumanRR.mat Chr1 '-v7.3' '-append'
   save WholeGenomeHumanRRE.mat counter '-append'
end
clear GM12878 SizeGM12878

%%
% This section is to call all the variable needed
% Cell line Embryonic stem cells

load('C:\Users\Biomecanica1\Documents\MATLAB\Final\SizeRR.mat', 'SizeH1esc') % Variable with the length of the RR genomic signals
load('C:\Users\Biomecanica1\Documents\MATLAB\DNA\H1esc.mat', 'H1esc') % Variable with the extracted RR DNA sequences

% In this section, the DNA sequence are mapped to a genomic signal and also the step one and two of the normalization
% is performed.

for l=1:15
   for m=1:size(H1esc(1).Region(l).Zone,2)

       R=find(upper(H1esc(1).Region(l).Zone(m).Sequence)=='N'); % To evaluate if a sequence has an any nucleotide within
       counter(1,l,2)=0;
           
       if (isempty(R))
               
           Signal=mapping(H1esc(1).Region(l).Zone(m).Sequence,3); % Neighbor joining mapping
           Fourier=abs(fft(Signal));
           Chr1(2,l).Matrix(m-counter(1,l,2),:)=interp1(1:size(Fourier,2),Fourier,linspace(1,size(Fourier,2),WholeMedian)); % Step one of the normalization
           Chr1(2,l).Matrix(m-counter(1,l,2),:)=Chr1(2,l).Matrix(m-counter(1,l,2),:)/SizeH1esc{1,l}(m); % Step two of the normalization
               
       else
               
           counter(1,l,2)=counter(1,l,2)+1; % To count if a sequence has an any nucleotide within
           
       end
           
       clear Signal Fourier R

   end
   display(l)
   save WholeGenomeHumanRR.mat Chr1 '-v7.3' '-append'
   save WholeGenomeHumanRRE.mat counter '-append'
end
clear H1esc SizeH1esc

%%
% This section is to call all the variable needed
% Cell line Hepatocellular carcinoma cells

load('C:\Users\Biomecanica1\Documents\MATLAB\Final\SizeRR.mat', 'SizeHepg2') % Variable with the length of the RR genomic signals
load('C:\Users\Biomecanica1\Documents\MATLAB\DNA\Hepg2.mat', 'Hepg2') % Variable with the extracted RR DNA sequences

% In this section, the DNA sequence are mapped to a genomic signal and also the step one and two of the normalization
% is performed.

for l=1:15
   for m=1:size(Hepg2(1).Region(l).Zone,2)

       R=find(upper(Hepg2(1).Region(l).Zone(m).Sequence)=='N'); % To evaluate if a sequence has an any nucleotide within
       counter(1,l,3)=0;
           
       if (isempty(R))
               
           Signal=mapping(Hepg2(1).Region(l).Zone(m).Sequence,3); % Neighbor joining mapping
           Fourier=abs(fft(Signal));
           Chr1(3,l).Matrix(m-counter(1,l,3),:)=interp1(1:size(Fourier,2),Fourier,linspace(1,size(Fourier,2),WholeMedian)); % Step one of the normalization
           Chr1(3,l).Matrix(m-counter(1,l,3),:)=Chr1(3,l).Matrix(m-counter(1,l,3),:)/SizeHepg2{1,l}(m); % Step two of the normalization
               
       else
               
           counter(1,l,3)=counter(1,l,3)+1; % To count if a sequence has an any nucleotide within
           
       end
           
       clear Signal Fourier R

   end
   display(l)
   save WholeGenomeHumanRR.mat Chr1 '-v7.3' '-append'
   save WholeGenomeHumanRRE.mat counter '-append'
end
clear Hepg2 SizeHepg2

%%
% This section is to call all the variable needed
% Cell line Mammary epithelial cells

load('C:\Users\Biomecanica1\Documents\MATLAB\Final\SizeRR.mat', 'SizeHmec') % Variable with the length of the RR genomic signals
load('C:\Users\Biomecanica1\Documents\MATLAB\DNA\Hmec.mat', 'Hmec') % Variable with the extracted RR DNA sequences

% In this section, the DNA sequence are mapped to a genomic signal and also the step one and two of the normalization
% is performed.

for l=1:15
   for m=1:size(Hmec(1).Region(l).Zone,2)

       R=find(upper(Hmec(1).Region(l).Zone(m).Sequence)=='N'); % To evaluate if a sequence has an any nucleotide within
       counter(1,l,4)=0;
           
       if (isempty(R))
               
           Signal=mapping(Hmec(1).Region(l).Zone(m).Sequence,3); % Neighbor joining mapping
           Fourier=abs(fft(Signal));
           Chr1(4,l).Matrix(m-counter(1,l,4),:)=interp1(1:size(Fourier,2),Fourier,linspace(1,size(Fourier,2),WholeMedian)); % Step one of the normalization
           Chr1(4,l).Matrix(m-counter(1,l,4),:)=Chr1(4,l).Matrix(m-counter(1,l,4),:)/SizeHmec{1,l}(m); % Step two of the normalization
               
       else
               
           counter(1,l,4)=counter(1,l,4)+1; % To count if a sequence has an any nucleotide within
            
       end
           
       clear Signal Fourier R

   end
   display(l)
   save WholeGenomeHumanRR.mat Chr1 '-v7.3' '-append'
   save WholeGenomeHumanRRE.mat counter '-append'
end
clear Hmec SizeHmec

%%
% This section is to call all the variable needed
% Cell line Skeletal muscle myoblasts

load('C:\Users\Biomecanica1\Documents\MATLAB\Final\SizeRR.mat', 'SizeHsmm') % Variable with the length of the RR genomic signals
load('C:\Users\Biomecanica1\Documents\MATLAB\DNA\Hsmm.mat', 'Hsmm') % Variable with the extracted RR DNA sequences

% In this section, the DNA sequence are mapped to a genomic signal and also the step one and two of the normalization
% is performed.

for l=1:15
   for m=1:size(Hsmm(1).Region(l).Zone,2)

       R=find(upper(Hsmm(1).Region(l).Zone(m).Sequence)=='N'); % To evaluate if a sequence has an any nucleotide within
       counter(1,l,5)=0;
           
       if (isempty(R))
               
           Signal=mapping(Hsmm(1).Region(l).Zone(m).Sequence,3); % Neighbor joining mapping
           Fourier=abs(fft(Signal));
           Chr1(5,l).Matrix(m-counter(1,l,5),:)=interp1(1:size(Fourier,2),Fourier,linspace(1,size(Fourier,2),WholeMedian)); % Step one of the normalization
           Chr1(5,l).Matrix(m-counter(1,l,5),:)=Chr1(5,l).Matrix(m-counter(1,l,5),:)/SizeHsmm{1,l}(m); % Step two of the normalization
               
       else
               
           counter(1,l,5)=counter(1,l,5)+1; % To count if a sequence has an any nucleotide within
           
       end
           
       clear Signal Fourier R

   end
   display(l)
   save WholeGenomeHumanRR.mat Chr1 '-v7.3' '-append'
   save WholeGenomeHumanRRE.mat counter '-append'
end
clear Hsmm SizeHsmm

%%
% This section is to call all the variable needed
% Cell line Umbilical vein endothelial cells

load('C:\Users\Biomecanica1\Documents\MATLAB\Final\SizeRR.mat', 'SizeHuevec') % Variable with the length of the RR genomic signals
load('C:\Users\Biomecanica1\Documents\MATLAB\DNA\Huevec.mat', 'Huevec') % Variable with the extracted RR DNA sequences

% In this section, the DNA sequence are mapped to a genomic signal and also the step one and two of the normalization
% is performed.

for l=1:15
   for m=1:size(Huevec(1).Region(l).Zone,2)

       R=find(upper(Huevec(1).Region(l).Zone(m).Sequence)=='N'); % To evaluate if a sequence has an any nucleotide within
       counter(1,l,6)=0;
           
       if (isempty(R))
               
           Signal=mapping(Huevec(1).Region(l).Zone(m).Sequence,3); % Neighbor joining mapping
           Fourier=abs(fft(Signal));
           Chr1(6,l).Matrix(m-counter(1,l,6),:)=interp1(1:size(Fourier,2),Fourier,linspace(1,size(Fourier,2),WholeMedian)); % Step one of the normalization
           Chr1(6,l).Matrix(m-counter(1,l,6),:)=Chr1(6,l).Matrix(m-counter(1,l,6),:)/SizeHuevec{1,l}(m); % Step two of the normalization
               
       else
               
           counter(1,l,6)=counter(1,l,6)+1; % To count if a sequence has an any nucleotide within
           
       end
           
       clear Signal Fourier R

   end
   display(l)
   save WholeGenomeHumanRR.mat Chr1 '-v7.3' '-append'
   save WholeGenomeHumanRRE.mat counter '-append'
end

clear Huevec SizeHuevec

%%
% This section is to call all the variable needed
% Cell line Erythrocytic leukemia cells

load('C:\Users\Biomecanica1\Documents\MATLAB\Final\SizeRR.mat', 'SizeK562') % Variable with the length of the RR genomic signals
load('C:\Users\Biomecanica1\Documents\MATLAB\DNA\K562.mat', 'K562') % Variable with the extracted RR DNA sequences

% In this section, the DNA sequence are mapped to a genomic signal and also the step one and two of the normalization
% is performed.

for l=1:15
   for m=1:size(K562(1).Region(l).Zone,2)

       R=find(upper(K562(1).Region(l).Zone(m).Sequence)=='N'); % To evaluate if a sequence has an any nucleotide within
       counter(1,l,7)=0;
           
       if (isempty(R))
               
           Signal=mapping(K562(1).Region(l).Zone(m).Sequence,3); % Neighbor joining mapping
           Fourier=abs(fft(Signal));
           Chr1(7,l).Matrix(m-counter(1,l,7),:)=interp1(1:size(Fourier,2),Fourier,linspace(1,size(Fourier,2),WholeMedian)); % Step one of the normalization
           Chr1(7,l).Matrix(m-counter(1,l,7),:)=Chr1(7,l).Matrix(m-counter(1,l,7),:)/SizeK562{1,l}(m); % Step two of the normalization
               
       else
               
           counter(1,l,7)=counter(1,l,7)+1; % To count if a sequence has an any nucleotide within
           
       end
           
       clear Signal Fourier R

   end
   display(l)
   save WholeGenomeHumanRR.mat Chr1 '-v7.3' '-append'
   save WholeGenomeHumanRRE.mat counter '-append'
end

clear K562 SizeK562

%%
% This section is to call all the variable needed
% Cell line Normal epidermal keratinocytes

load('C:\Users\Biomecanica1\Documents\MATLAB\Final\SizeRR.mat', 'SizeNhek') % Variable with the length of the RR genomic signals
load('C:\Users\Biomecanica1\Documents\MATLAB\DNA\Nhek.mat', 'Nhek') % Variable with the extracted RR DNA sequences

% In this section, the DNA sequence are mapped to a genomic signal and also the step one and two of the normalization
% is performed.

for l=1:15
   for m=1:size(Nhek(1).Region(l).Zone,2)

       R=find(upper(Nhek(1).Region(l).Zone(m).Sequence)=='N'); % To evaluate if a sequence has an any nucleotide within
       counter(1,l,8)=0;
           
       if (isempty(R))
               
           Signal=mapping(Nhek(1).Region(l).Zone(m).Sequence,3); % Neighbor joining mapping
           Fourier=abs(fft(Signal));
           Chr1(8,l).Matrix(m-counter(1,l,8),:)=interp1(1:size(Fourier,2),Fourier,linspace(1,size(Fourier,2),WholeMedian)); % Step one of the normalization
           Chr1(8,l).Matrix(m-counter(1,l,8),:)=Chr1(8,l).Matrix(m-counter(1,l,8),:)/SizeNhek{1,l}(m); % Step two of the normalization
               
       else
               
           counter(1,l,8)=counter(1,l,8)+1; % To count if a sequence has an any nucleotide within
           
       end
           
       clear Signal Fourier R

   end
   display(l)
   save WholeGenomeHumanRR.mat Chr1 '-v7.3' '-append'
   save WholeGenomeHumanRRE.mat counter '-append'
end

clear Nhek SizeNhek

%%
% This section is to call all the variable needed
% Cell line Normal lung fibroblasts

load('C:\Users\Biomecanica1\Documents\MATLAB\Final\SizeRR.mat', 'SizeNhlf') % Variable with the length of the RR genomic signals
load('C:\Users\Biomecanica1\Documents\MATLAB\DNA\Nhlf.mat', 'Nhlf') % Variable with the extracted RR DNA sequences

% In this section, the DNA sequence are mapped to a genomic signal and also the step one and two of the normalization
% is performed.

for l=1:15
   for m=1:size(Nhlf(1).Region(l).Zone,2)

       R=find(upper(Nhlf(1).Region(l).Zone(m).Sequence)=='N'); % To evaluate if a sequence has an any nucleotide within
       counter(1,l,9)=0;
           
       if (isempty(R))
               
           Signal=mapping(Nhlf(1).Region(l).Zone(m).Sequence,3); % Neighbor joining mapping
           Fourier=abs(fft(Signal));
           Chr1(9,l).Matrix(m-counter(1,l,9),:)=interp1(1:size(Fourier,2),Fourier,linspace(1,size(Fourier,2),WholeMedian)); % Step one of the normalization
           Chr1(9,l).Matrix(m-counter(1,l,9),:)=Chr1(9,l).Matrix(m-counter(1,l,9),:)/SizeNhlf{1,l}(m); % Step two of the normalization
                
       else
               
           counter(1,l,9)=counter(1,l,9)+1; % To count if a sequence has an any nucleotide within
           
       end
           
       clear Signal Fourier R

   end
   display(l)
   save WholeGenomeHumanRR.mat Chr1 '-v7.3' '-append'
   save WholeGenomeHumanRRE.mat counter '-append'
end

clear Nhlf SizeNhlf
