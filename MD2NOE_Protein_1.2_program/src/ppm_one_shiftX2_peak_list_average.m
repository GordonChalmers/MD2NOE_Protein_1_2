

total_peak_list=cell(size(peak_list_shiftx,1),size(peak_list_shiftx,2));

for i=1:6:6*size(residues,2)
    
    maximum_ppm_one=
    
    for j=1:size(peak_list_shiftx(:,i),1)

        if isnan(peak_list_shiftx{j,i+2})==1 
            peak_list_shiftx{j,i+2}=999;
        end
        
        if isnan(peak_list_ppm_one{j,i})==1
            peak_list_ppm_one{j,i+2}=999;
        end
        
        if peak_list_shiftx{j,i+2}~=999
            if peak_list_ppm_one{j,i+2}~=999
                total_peak_list{j,i}=peak_list_shiftx{j,i};
                total_peak_list{j,i+1}=peak_list_shiftx{j,i+1};
                total_peak_list{j,i+2}=(peak_list_ppm_one{j,i+2}+peak_list_shiftx{j,i+2})/2;
                total_peak_list{j,i+3}=peak_list_shiftx{j,i+3};
                total_peak_list{j,i+4}=peak_list_shiftx{j,i+4};
                total_peak_list{j,i+5}=peak_list_shiftx{j,i+5};
            end
        end
        
        if peak_list_shiftx{j,i+2}~=999
            if peak_list_ppm_one{j,i+2}==999
                total_peak_list{j,i}=peak_list_shiftx{j,i};
                total_peak_list{j,i+1}=peak_list_shiftx{j,i+1};
                total_peak_list{j,i+2}=peak_list_shiftx{j,i+2};
                total_peak_list{j,i+3}=peak_list_shiftx{j,i+3};
                total_peak_list{j,i+4}=peak_list_shiftx{j,i+4};
                total_peak_list{j,i+5}=peak_list_shiftx{j,i+5};
            end
        end
        
        if peak_list_shiftx{j,i+2}==999
            if peak_list_ppm_one{j,i+2}~=999
                total_peak_list{j,i}=peak_list_ppm_one{j,i};
                total_peak_list{j,i+1}=peak_list_ppm_one{j,i+1};
                total_peak_list{j,i+2}=peak_list_ppm_one{j,i+2};
                total_peak_list{j,i+3}=peak_list_ppm_one{j,i+3};
                total_peak_list{j,i+4}=peak_list_ppm_one{j,i+4};
                total_peak_list{j,i+5}=peak_list_ppm_one{j,i+5};
            end
        end
        
        if peak_list_shiftx{j,i+2}==999
            if peak_list_ppm_one{j,i+2}==999
                total_peak_list{j,i}=peak_list_shiftx{j,i};
                total_peak_list{j,i+1}=peak_list_shiftx{j,i+1};
                total_peak_list{j,i+2}=(peak_list_ppm_one{j,i+2}+peak_list_shiftx{j,i+2})/2;
                total_peak_list{j,i+3}=peak_list_shiftx{j,i+3};
                total_peak_list{j,i+4}=peak_list_shiftx{j,i+4};
                total_peak_list{j,i+5}=peak_list_shiftx{j,i+5};
            end
        end
        
        
    end
    
end