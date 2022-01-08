
%% This program takes in a peak list file and generates
%%  an excel sheet of noe vectors.

%% Inputs are the peak list file, the width of the peaks,
%%  and the chemical shift list.

%% input_peak_list_file
%% chemical_shift_width
%% output_file_spectrum

global width;
global peak_list
global chemical_shift_list;

%% width=chemical_shift_width;
%% chemical_shift_list=csvread(input_chemical_shift_file);
%% file_spectrumID=fopen(output_file_spectrum,'w');

%% input peak list file 
%% set width, e.g. .2 ppm
%% chemical shift list is the set of chemical shifts in the experimental data 


max_peak=zeros(16,1);

%% residues=[132 155 157 171 208 216 240 274 275 290 340 356 357 371 390 398];]

%% input_peak_list_file="st6_2ndRun_calculations/1000ns/peak_list_1_1000_2ndRun_not_averaged.xlsx";

%% peak_list=xlsread(input_peak_list_file);

size_residues=size(peak_list,2)/6;

spectra=zeros(size(chemical_shift_list,1),size_residues);

for i=1:size(peak_list,2)/6
    
    for i2=1:size(peak_list,1)
    
        if peak_list{i2,6*(i-1)+3}~=999 

           if isnan(peak_list{i2,6*(i-1)+3})~=1
            
                for j=1:size(chemical_shift_list,1)

                    spectra(j,i)=spectra(j,i)+peak_list{i2,6*(i-1)+5}*exp(-(peak_list{i2,6*(i-1)+3}-chemical_shift_list(j))^2/2/width^2);

                end
           end

                
        end
        
    end
    
end


%% add autopeak

for i=1:(size(peak_list,2))/6

    max_peak(i)=min(spectra(:,i));
    
end
    
max_peak_spectrum=mean(max_peak(:));

for i=1:(size(peak_list,2))/6
    
    for j=1:size(chemical_shift_list,1)
        spectra(j,i)=spectra(j,i)+max_peak_spectrum*exp(-(chemical_shift_list(size(chemical_shift_list,1))-chemical_shift_list(j))^2/2/width^2);
    end
    
end


%% save the peak list and spectra

save('shiftx2_output_peak_list_spectra','peak_list','spectra');





