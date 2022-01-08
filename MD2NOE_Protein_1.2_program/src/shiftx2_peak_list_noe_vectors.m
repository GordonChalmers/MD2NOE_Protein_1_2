

%% the shiftX2 trajectory output file is in a different format than the 
%%  the single frame calculation.
%%
%% this script calculates the spectra from the shiftX2 chemical shift 
%%  file
%%
%% it has to be changed for different proteins
%%


max_peak=zeros(16,1);

residues=[29,52,54,68,105,113,137,171,172,187,237,253,254,268,287,295];

spectra=zeros(size(chemical_shift_list,1),16);



    for i2=1:size(peak_list,1)
        

if peak_list(i2,1)==29  
   i=1;
end

if peak_list(i2,1)==52 
   i=2;
end

if peak_list(i2,1)==54                 
   i=3;
end

if peak_list(i2,1)==68 
   i=4;
end

if peak_list(i2,1)==105  
   i=5;
end1

if peak_list(i2,1)==113 
   i=6;
end

if peak_list(i2,1)==137 
   i=7;
end

if peak_list(i2,1)==171 
   i=8;
end

if peak_list(i2,1)==172  
   i=9;
end1

if peak_list(i2,1)==187
   i=10;
end

if peak_list(i2,1)==237 
   i=11;
end

if peak_list(i2,1)==253 
   i=12;
end

if peak_list(i2,1)==254  
   i=13;
end

if peak_list(i2,1)==268 
   i=14;
end

if peak_list(i2,1)==287 
   i=15;
end

if peak_list(i2,1)==295 
   i=16;
end

              if peak_list(i2,2)~=999

                for j=1:size(chemical_shift_list,1)
                    spectra(j,i)=spectra(j,i)+peak_list(i2,3)*exp(-(peak_list(i2,2)-chemical_shift_list(j))^2/2/.2^2);
                end
                
                if peak_list(i2,3)<max_peak(i)
                    max_peak(i)=peak_list(i2,3);
                end

             end

                
    end
            

   
%% add autopeak

for i=1:size(peak_list,2)/6

    max_peak(i)=min(spectra(:,i));
    
end
    
max_peak_spectrum=mean(max_peak(:));

for i=1:size(peak_list,2)/6
    
    for j=1:size(chemical_shift_list,1)
        spectra(j,i)=spectra(j,i)+max_peak_spectrum*exp(-(chemical_shift_list(size(chemical_shift_list,1))-chemical_shift_list(j))^2/2/width^2);
    end
    
end

            
save('output_file_spectrum','peak_list','spectra');




