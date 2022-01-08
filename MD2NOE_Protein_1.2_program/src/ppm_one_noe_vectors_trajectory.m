
%% Inputs are a ppm_one file of chemical shifts

%% The output is a peak list which can be used in noe_vectors_from_peak_list.m.

%% ppm1 chemical shift files

chemical_shift_ppm_one_fileID=fopen(input_chemical_ppm_one_file,'r');
chemical_shift_ppm_one_second_fileID=fopen(input_chemical_ppm_one_second_file,'r');

%% noe information file

noe_input_fileID=fopen(noe_input_file,'r');


mixing=mixing_milliseconds;

residues=input_residues;


%% at most 100 peaks max_distance

peak_list=cell(10,6*size(residues,2));


%% residues from input

for i=1:size(residues,2)
    
    total=0;
    frewind(noe_input_fileID)
    
    %% check for all H- spin pairs in noe's
    
    while ~feof(noe_input_fileID)
        noe_file=fgetl(noe_input_fileID);
        if noe_file~=-2
            test=strsplit(noe_file,',');
            
            if str2double(test{1,3})==residues(i)
                if strcmp(test{1,4},'H')
                    
                    %% atom in spin pair from residue i
                    
                    spin_pair_atom=test{1,7};
                    
                    spin_pair_residue_type=test{1,5};
                    spin_pair_residue=str2double(test{1,6});
                    
                    total=total+1;
                    
                    peak_list{total,6*(i-1)+1}=spin_pair_residue;
                    peak_list{total,6*(i-1)+2}=spin_pair_atom;
                    peak_list{total,6*(i-1)+4}=spin_pair_residue_type;
                    
                    %% parse the noe files
                    
                    filename=test{1,10};
                    
                    noe_fileID=fopen(filename,'r');
                    
                    while ~feof(noe_fileID)
                        
                        noe_test=fgetl(noe_fileID);
                        
                        noe_test_split=strsplit(noe_test,' ');
                        
                        if str2double(noe_test_split{1,1})==mixing
                            
                            noe=str2double(noe_test_split{1,2});
                            
                        end
                        
                    end
                    
                    fclose(noe_fileID);
                    
                    peak_list{total,6*(i-1)+5}=noe;
                    
                    
                    peak_list{total,6*(i-1)+3}=999;
                    
                    %% find chemical shift
                    
                    
                    frewind(chemical_shift_ppm_one_fileID);
                    chem_file=fgetl(chemical_shift_ppm_one_fileID);
                    
                    no_match=1;
                    
                    while ~feof(chemical_shift_ppm_one_fileID)
                        chem_file=fgetl(chemical_shift_ppm_one_fileID);
                        if chem_file~=-1
                            test=strsplit(chem_file,' ');
                            
                            if strcmp(spin_pair_atom,test{1,4})==1
                                
                                if spin_pair_residue==str2double(test{1,2})
                                    
                                    peak_list{total,6*(i-1)+3}=str2double(test{1,5});
                                    
                                end
                                
                            end
                            
                            
                            if strcmp(spin_pair_atom,'HD1')==1
                                
                                if strcmp(spin_pair_residue_type,"ALA")==0&strcmp(spin_pair_residue_type,"THR")==0 ...
                                        &strcmp(spin_pair_residue_type,"LEU")==0&strcmp(spin_pair_residue_type,"ILE")==0 ...
                                        &strcmp(spin_pair_residue_type,"MET")==0&strcmp(spin_pair_residue_type,"VAL")==0
                                    
                                    if strcmp(test{1,4},'HD1')==1
                                        
                                        if spin_pair_residue==str2double(test{1,2})
                                            
                                            peak_list{total,6*(i-1)+3}=str2double(test{1,5});
                                            
                                        end
                                        
                                    end
                                    
                                end
                                
                            end
                            
                            
                            if strcmp(spin_pair_atom,'HD2')==1
                                
                                if strcmp(spin_pair_residue_type,"ALA")==0&strcmp(spin_pair_residue_type,"THR")==0 ...
                                        &strcmp(spin_pair_residue_type,"LEU")==0&strcmp(spin_pair_residue_type,"ILE")==0 ...
                                        &strcmp(spin_pair_residue_type,"MET")==0&strcmp(spin_pair_residue_type,"VAL")==0
                                    
                                    if strcmp(test{1,4},'HD2')==1
                                        
                                        if spin_pair_residue==str2double(test{1,2})
                                            
                                            peak_list{total,6*(i-1)+3}=str2double(test{1,5});
                                            
                                        end
                                        
                                    end
                                    
                                end
                                
                            end
                            
                            if strcmp(spin_pair_atom,'HD3')==1
                                
                                if strcmp(spin_pair_residue_type,"ALA")==0&strcmp(spin_pair_residue_type,"THR")==0 ...
                                        &strcmp(spin_pair_residue_type,"LEU")==0&strcmp(spin_pair_residue_type,"ILE")==0 ...
                                        &strcmp(spin_pair_residue_type,"MET")==0&strcmp(spin_pair_residue_type,"VAL")==0
                                    
                                    if strcmp(test{1,4},'HD3')==1
                                        
                                        if spin_pair_residue==str2double(test{1,2})
                                            
                                            peak_list{total,6*(i-1)+3}=str2double(test{1,5});
                                            
                                        end
                                        
                                    end
                                    
                                end
                                
                            end
                            
                            if strcmp(spin_pair_atom,'HD1')==1
                                
                                if strcmp(spin_pair_residue_type,"ALA")==0&strcmp(spin_pair_residue_type,"THR")==0 ...
                                        &strcmp(spin_pair_residue_type,"LEU")==0&strcmp(spin_pair_residue_type,"ILE")==0 ...
                                        &strcmp(spin_pair_residue_type,"MET")==0&strcmp(spin_pair_residue_type,"VAL")==0
                                    
                                    if strcmp(test{1,4},'HD')==1
                                        
                                        if spin_pair_residue==str2double(test{1,2})
                                            
                                            peak_list{total,6*(i-1)+3}=str2double(test{1,5});
                                            
                                        end
                                        
                                    end
                                    
                                end
                                
                            end
                            
                            
                            if strcmp(spin_pair_atom,'HD2')==1
                                
                                if strcmp(spin_pair_residue_type,"ALA")==0&strcmp(spin_pair_residue_type,"THR")==0 ...
                                        &strcmp(spin_pair_residue_type,"LEU")==0&strcmp(spin_pair_residue_type,"ILE")==0 ...
                                        &strcmp(spin_pair_residue_type,"MET")==0&strcmp(spin_pair_residue_type,"VAL")==0
                                    
                                    if strcmp(test{1,4},'HD')==1
                                        
                                        if spin_pair_residue==str2double(test{1,2})
                                            
                                            peak_list{total,6*(i-1)+3}=str2double(test{1,5});
                                            
                                        end
                                        
                                    end
                                    
                                end
                                
                            end
                            
                            
                            if strcmp(spin_pair_atom,'HG2')==1
                                
                                if strcmp(spin_pair_residue_type,"ALA")==0&strcmp(spin_pair_residue_type,"THR")==0 ...
                                        &strcmp(spin_pair_residue_type,"LEU")==0&strcmp(spin_pair_residue_type,"ILE")==0 ...
                                        &strcmp(spin_pair_residue_type,"MET")==0&strcmp(spin_pair_residue_type,"VAL")==0
                                    
                                    if strcmp(test{1,4},'HG')==1
                                        
                                        if spin_pair_residue==str2double(test{1,2})
                                            
                                            peak_list{total,6*(i-1)+3}=str2double(test{1,5});
                                            
                                        end
                                        
                                    end
                                    
                                end
                                
                            end
                            
                            
                            if strcmp(spin_pair_atom,'HG3')==1
                                
                                if strcmp(spin_pair_residue_type,"ALA")==0&strcmp(spin_pair_residue_type,"THR")==0 ...
                                        &strcmp(spin_pair_residue_type,"LEU")==0&strcmp(spin_pair_residue_type,"ILE")==0 ...
                                        &strcmp(spin_pair_residue_type,"MET")==0&strcmp(spin_pair_residue_type,"VAL")==0
                                    
                                    if strcmp(test{1,4},'HG')==1
                                        
                                        if spin_pair_residue==str2double(test{1,2})
                                            
                                            peak_list{total,6*(i-1)+3}=str2double(test{1,5});
                                            
                                        end
                                        
                                    end
                                    
                                end
                                
                            end
                            
                            
                            
                            if strcmp(spin_pair_atom,'HE1')==1
                                
                                if strcmp(spin_pair_residue_type,"ALA")==0&strcmp(spin_pair_residue_type,"THR")==0 ...
                                        &strcmp(spin_pair_residue_type,"LEU")==0&strcmp(spin_pair_residue_type,"ILE")==0 ...
                                        &strcmp(spin_pair_residue_type,"MET")==0&strcmp(spin_pair_residue_type,"VAL")==0
                                    
                                    if strcmp(test{1,4},'HE')==1
                                        
                                        if spin_pair_residue==str2double(test{1,2})
                                            
                                            peak_list{total,6*(i-1)+3}=str2double(test{1,5});
                                            
                                        end
                                        
                                    end
                                    
                                end
                                
                            end
                            
                            
                            if strcmp(spin_pair_atom,'HE2')==1
                                
                                if strcmp(spin_pair_residue_type,"ALA")==0&strcmp(spin_pair_residue_type,"THR")==0 ...
                                        &strcmp(spin_pair_residue_type,"LEU")==0&strcmp(spin_pair_residue_type,"ILE")==0 ...
                                        &strcmp(spin_pair_residue_type,"MET")==0&strcmp(spin_pair_residue_type,"VAL")==0
                                    
                                    if strcmp(test{1,4},'HE')==1
                                        
                                        if spin_pair_residue==str2double(test{1,2})
                                            
                                            peak_list{total,6*(i-1)+3}=str2double(test{1,5});
                                            
                                        end
                                        
                                    end
                                    
                                end
                                
                            end
                            
                            
                            
                            if strcmp(spin_pair_atom,'HG11')==1
                                
                                if strcmp(spin_pair_residue_type,"ALA")==1|strcmp(spin_pair_residue_type,"THR")==1 ...
                                        |strcmp(spin_pair_residue_type,"LEU")==1|strcmp(spin_pair_residue_type,"ILE")==1 ...
                                        |strcmp(spin_pair_residue_type,"MET")==1|strcmp(spin_pair_residue_type,"VAL")==1
                                    
                                    if strcmp(test{1,4},'HG1')==1
                                        
                                        if spin_pair_residue==str2double(test{1,2})
                                            
                                            peak_list{total,6*(i-1)+3}=str2double(test{1,5});
                                            
                                        end
                                        
                                    end
                                    
                                end
                                
                            end
                            
                            if strcmp(spin_pair_atom,'HG12')==1
                                
                                if strcmp(spin_pair_residue_type,"ALA")==1|strcmp(spin_pair_residue_type,"THR")==1 ...
                                        |strcmp(spin_pair_residue_type,"LEU")==1|strcmp(spin_pair_residue_type,"ILE")==1 ...
                                        |strcmp(spin_pair_residue_type,"MET")==1|strcmp(spin_pair_residue_type,"VAL")==1
                                    
                                    if strcmp(test{1,4},'HG1')==1
                                        
                                        if spin_pair_residue==str2double(test{1,2})
                                            
                                            peak_list{total,6*(i-1)+3}=str2double(test{1,5});
                                            
                                        end
                                        
                                    end
                                    
                                end
                                
                            end
                            
                            
                            if strcmp(spin_pair_atom,'HG13')==1
                                
                                if strcmp(spin_pair_residue_type,"ALA")==1|strcmp(spin_pair_residue_type,"THR")==1 ...
                                        |strcmp(spin_pair_residue_type,"LEU")==1|strcmp(spin_pair_residue_type,"ILE")==1 ...
                                        |strcmp(spin_pair_residue_type,"MET")==1|strcmp(spin_pair_residue_type,"VAL")==1
                                    
                                    if strcmp(test{1,4},'HG1')==1
                                        
                                        if spin_pair_residue==str2double(test{1,2})
                                            
                                            peak_list{total,6*(i-1)+3}=str2double(test{1,5});
                                            
                                        end
                                        
                                    end
                                    
                                end
                                
                            end
                            
                            
                            if strcmp(spin_pair_atom,'HG21')==1
                                
                                if strcmp(spin_pair_residue_type,"ALA")==1|strcmp(spin_pair_residue_type,"THR")==1 ...
                                        |strcmp(spin_pair_residue_type,"LEU")==1|strcmp(spin_pair_residue_type,"ILE")==1 ...
                                        |strcmp(spin_pair_residue_type,"MET")==1|strcmp(spin_pair_residue_type,"VAL")==1
                                    
                                    
                                    if strcmp(test{1,4},'HG2')==1
                                        
                                        if spin_pair_residue==str2double(test{1,2})
                                            
                                            peak_list{total,6*(i-1)+3}=str2double(test{1,5});
                                            
                                        end
                                        
                                    end
                                    
                                end
                                
                            end
                            
                            
                            if strcmp(spin_pair_atom,'HG22')==1
                                
                                if strcmp(spin_pair_residue_type,"ALA")==1|strcmp(spin_pair_residue_type,"THR")==1 ...
                                        |strcmp(spin_pair_residue_type,"LEU")==1|strcmp(spin_pair_residue_type,"ILE")==1 ...
                                        |strcmp(spin_pair_residue_type,"MET")==1|strcmp(spin_pair_residue_type,"VAL")==1
                                    
                                    if strcmp(test{1,4},'HG2')==1
                                        
                                        if spin_pair_residue==str2double(test{1,2})
                                            
                                            peak_list{total,6*(i-1)+3}=str2double(test{1,5});
                                            
                                        end
                                        
                                    end
                                    
                                end
                                
                            end
                            
                            
                            if strcmp(spin_pair_atom,'HG23')==1
                                
                                if strcmp(spin_pair_residue_type,"ALA")==1|strcmp(spin_pair_residue_type,"THR")==1 ...
                                        |strcmp(spin_pair_residue_type,"LEU")==1|strcmp(spin_pair_residue_type,"ILE")==1 ...
                                        |strcmp(spin_pair_residue_type,"MET")==1|strcmp(spin_pair_residue_type,"VAL")==1
                                    
                                    if strcmp(test{1,4},'HG2')==1
                                        
                                        if spin_pair_residue==str2double(test{1,2})
                                            
                                            peak_list{total,6*(i-1)+3}=str2double(test{1,5});
                                            
                                        end
                                        
                                    end
                                    
                                end
                                
                            end
                            
                            
                            if strcmp(spin_pair_atom,'HD11')==1
                                
                                if strcmp(spin_pair_residue_type,"ALA")==1|strcmp(spin_pair_residue_type,"THR")==1 ...
                                        |strcmp(spin_pair_residue_type,"LEU")==1|strcmp(spin_pair_residue_type,"ILE")==1 ...
                                        |strcmp(spin_pair_residue_type,"MET")==1|strcmp(spin_pair_residue_type,"VAL")==1
                                    
                                    if strcmp(test{1,4},'HD1')==1
                                        
                                        if spin_pair_residue==str2double(test{1,2})
                                            
                                            peak_list{total,6*(i-1)+3}=str2double(test{1,5});
                                            
                                        end
                                        
                                    end
                                    
                                end
                                
                            end
                            
                            if strcmp(spin_pair_atom,'HD12')==1
                                
                                if strcmp(spin_pair_residue_type,"ALA")==1|strcmp(spin_pair_residue_type,"THR")==1 ...
                                        |strcmp(spin_pair_residue_type,"LEU")==1|strcmp(spin_pair_residue_type,"ILE")==1 ...
                                        |strcmp(spin_pair_residue_type,"MET")==1|strcmp(spin_pair_residue_type,"VAL")==1
                                    
                                    if strcmp(test{1,4},'HD1')==1
                                        
                                        if spin_pair_residue==str2double(test{1,2})
                                            
                                            peak_list{total,6*(i-1)+3}=str2double(test{1,5});
                                            
                                        end
                                        
                                    end
                                    
                                end
                                
                            end
                            
                            
                            if strcmp(spin_pair_atom,'HD13')==1
                                
                                if strcmp(spin_pair_residue_type,"ALA")==1|strcmp(spin_pair_residue_type,"THR")==1 ...
                                        |strcmp(spin_pair_residue_type,"LEU")==1|strcmp(spin_pair_residue_type,"ILE")==1 ...
                                        |strcmp(spin_pair_residue_type,"MET")==1|strcmp(spin_pair_residue_type,"VAL")==1
                                    
                                    if strcmp(test{1,4},'HD1')==1
                                        
                                        if spin_pair_residue==str2double(test{1,2})
                                            
                                            peak_list{total,6*(i-1)+3}=str2double(test{1,5});
                                            
                                        end
                                        
                                    end
                                    
                                end
                                
                            end
                            
                            
                            if strcmp(spin_pair_atom,'HE1')==1
                                
                                if strcmp(spin_pair_residue_type,"MET")==1
                                    
                                    if strcmp(test{1,4},'HE')==1
                                        
                                        if spin_pair_residue==str2double(test{1,2})
                                            
                                            peak_list{total,6*(i-1)+3}=str2double(test{1,5});
                                            
                                        end
                                        
                                    end
                                    
                                end
                                
                            end
                            
                            if strcmp(spin_pair_atom,'HE2')==1
                                
                                if strcmp(spin_pair_residue_type,"MET")==1
                                    
                                    if strcmp(test{1,4},'HE')==1
                                        
                                        if spin_pair_residue==str2double(test{1,2})
                                            
                                            peak_list{total,6*(i-1)+3}=str2double(test{1,5});
                                            
                                        end
                                        
                                    end
                                    
                                end
                                
                            end
                            
                            
                            if strcmp(spin_pair_atom,'HE3')==1
                                
                                if strcmp(spin_pair_residue_type,"MET")==1
                                    
                                    if strcmp(test{1,4},'HE')==1
                                        
                                        if spin_pair_residue==str2double(test{1,2})
                                            
                                            peak_list{total,6*(i-1)+3}=str2double(test{1,5});
                                            
                                        end
                                        
                                    end
                                    
                                end
                                
                            end
                            
                            
                            if strcmp(spin_pair_atom,'HB1')==1
                                
                                if strcmp(spin_pair_residue_type,"ALA")==1
                                    
                                    if strcmp(test{1,4},'HB')==1
                                        
                                        if spin_pair_residue==str2double(test{1,2})
                                            
                                            peak_list{total,6*(i-1)+3}=str2double(test{1,5});
                                            
                                        end
                                        
                                    end
                                    
                                end
                                
                            end
                            
                            if strcmp(spin_pair_atom,'HB2')==1
                                
                                if strcmp(spin_pair_residue_type,"ALA")==1
                                    
                                    if strcmp(test{1,4},'HB')==1
                                        
                                        if spin_pair_residue==str2double(test{1,2})
                                            
                                            peak_list{total,6*(i-1)+3}=str2double(test{1,5});
                                            
                                        end
                                        
                                    end
                                    
                                end
                                
                            end
                            
                            
                            if strcmp(spin_pair_atom,'HB3')==1
                                
                                if strcmp(spin_pair_residue_type,"ALA")==1|strcmp(spin_pair_residue_type,"THR")==1 ...
                                        |strcmp(spin_pair_residue_type,"LEU")==1|strcmp(spin_pair_residue_type,"ILE")==1 ...
                                        |strcmp(spin_pair_residue_type,"MET")==1|strcmp(spin_pair_residue_type,"VAL")==1
                                    
                                    if strcmp(test{1,4},'HB')==1
                                        
                                        if spin_pair_residue==str2double(test{1,2})
                                            
                                            peak_list{total,6*(i-1)+3}=str2double(test{1,5});
                                            
                                        end
                                        
                                    end
                                    
                                end
                                
                            end
                            
                            
                            
                            if strcmp(spin_pair_atom,'HD21')==1
                                
                                if strcmp(spin_pair_residue_type,"ALA")==1|strcmp(spin_pair_residue_type,"THR")==1 ...
                                        |strcmp(spin_pair_residue_type,"LEU")==1|strcmp(spin_pair_residue_type,"ILE")==1 ...
                                        |strcmp(spin_pair_residue_type,"MET")==1|strcmp(spin_pair_residue_type,"VAL")==1
                                    
                                    
                                    if strcmp(test{1,4},'HD2')==1
                                        
                                        if spin_pair_residue==str2double(test{1,2})
                                            
                                            peak_list{total,6*(i-1)+3}=str2double(test{1,5});
                                            
                                        end
                                        
                                    end
                                    
                                end
                                
                            end
                            
                            
                            if strcmp(spin_pair_atom,'HD22')==1
                                
                                if strcmp(spin_pair_residue_type,"ALA")==1|strcmp(spin_pair_residue_type,"THR")==1 ...
                                        |strcmp(spin_pair_residue_type,"LEU")==1|strcmp(spin_pair_residue_type,"ILE")==1 ...
                                        |strcmp(spin_pair_residue_type,"MET")==1|strcmp(spin_pair_residue_type,"VAL")==1
                                    
                                    if strcmp(test{1,4},'HD2')==1
                                        
                                        if spin_pair_residue==str2double(test{1,2})
                                            
                                            peak_list{total,6*(i-1)+3}=str2double(test{1,5});
                                            
                                        end
                                        
                                    end
                                    
                                end
                                
                            end
                            
                            
                            if strcmp(spin_pair_atom,'HD23')==1
                                
                                if strcmp(spin_pair_residue_type,"ALA")==1|strcmp(spin_pair_residue_type,"THR")==1 ...
                                        |strcmp(spin_pair_residue_type,"LEU")==1|strcmp(spin_pair_residue_type,"ILE")==1 ...
                                        |strcmp(spin_pair_residue_type,"MET")==1|strcmp(spin_pair_residue_type,"VAL")==1
                                    
                                    if strcmp(test{1,4},'HD2')==1
                                        
                                        if spin_pair_residue==str2double(test{1,2})
                                            
                                            peak_list{total,6*(i-1)+3}=str2double(test{1,5});
                                            
                                        end
                                        
                                    end
                                    
                                end
                                
                            end
                            
                        end
                        
                    end
                    
                    total_test=0;
                    
                    frewind(chemical_shift_ppm_one_second_fileID);
                    chem_file=fgetl(chemical_shift_ppm_one_second_fileID);
                    
                    while ~feof(chemical_shift_ppm_one_second_fileID)
                        
                        chem_file=fgetl(chemical_shift_ppm_one_second_fileID);
                        if chem_file~=-1
                            test=strsplit(chem_file,{' ','-'});
                            
                            if strcmp(spin_pair_atom,'H')==1
                                if spin_pair_residue==str2double(test{1,2})
                                    
                                    peak_list{total,6*(i-1)+3}=str2double(test{1,10});
                                    
                                end
                                
                            end
                            
                            if strcmp(spin_pair_atom,'HA')==1
                                if spin_pair_residue==str2double(test{1,2})
                                    
                                    peak_list{total,6*(i-1)+3}=str2double(test{1,14});
                                    
                                end
                            end
                            
                            
                            if strcmp(spin_pair_atom,'HH')==1
                                
                                if strcmp(spin_pair_residue_type,"ALA")==0&strcmp(spin_pair_residue_type,"THR")==0 ...
                                        &strcmp(spin_pair_residue_type,"LEU")==0&strcmp(spin_pair_residue_type,"ILE")==0 ...
                                        &strcmp(spin_pair_residue_type,"MET")==0&strcmp(spin_pair_residue_type,"VAL")==0
                                    
                                    if spin_pair_residue==str2double(test{1,2})
                                        
                                        peak_list{total,6*(i-1)+3}=str2double(test{1,10});
                                        
                                    end
                                    
                                end
                                
                            end
                            
                            if strcmp(spin_pair_atom,'HA2')==1
                                
                                if strcmp(spin_pair_residue_type,"ALA")==0&strcmp(spin_pair_residue_type,"THR")==0 ...
                                        &strcmp(spin_pair_residue_type,"LEU")==0&strcmp(spin_pair_residue_type,"ILE")==0 ...
                                        &strcmp(spin_pair_residue_type,"MET")==0&strcmp(spin_pair_residue_type,"VAL")==0
                                    
                                    if spin_pair_residue==str2double(test{1,2})
                                        
                                        peak_list{total,6*(i-1)+3}=str2double(test{1,14});
                                        
                                    end
                                    
                                end
                                
                            end
                            
                            
                            if strcmp(spin_pair_atom,'HA3')==1
                                
                                if strcmp(spin_pair_residue_type,"ALA")==0&strcmp(spin_pair_residue_type,"THR")==0 ...
                                        &strcmp(spin_pair_residue_type,"LEU")==0&strcmp(spin_pair_residue_type,"ILE")==0 ...
                                        &strcmp(spin_pair_residue_type,"MET")==0&strcmp(spin_pair_residue_type,"VAL")==0
                                    
                                    if spin_pair_residue==str2double(test{1,2})
                                        
                                        peak_list{total,6*(i-1)+3}=str2double(test{1,14});
                                        
                                    end
                                    
                                end
                                
                            end
                            
                        end
                                               
                        
                    end
                end
            end
            
        end
    end
end

fclose(chemical_shift_ppm_one_fileID);
fclose(noe_input_fileID);
fclose(chemical_shift_ppm_one_second_fileID);


