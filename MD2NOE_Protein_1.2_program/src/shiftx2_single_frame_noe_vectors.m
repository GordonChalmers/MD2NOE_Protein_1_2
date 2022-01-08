
%% This program creates an noe spectrum from the 1/r^6 calculation.

%% Inputs are a shiftx2 file of chemical shifts, 1/r^6
%%  calculations from the SingleFrameNMR program, and an output
%%  file.

%% The output is a peak list which can be used in noe_vectors_from_peak_list.m.

chemical_shiftx2_fileID=fopen(input_chemical_shiftx2_file,'r');
noe_input_fileID=fopen(input_noe_single_frame_file,'r');

residues=input_residues;

%% at most 100 peaks max_distance

peak_list=cell(10,6*size(residues,2));

%% residues from input

for i=1:size(residues,2)
    
    total=0;
    frewind(noe_input_fileID)
    
    %% check for all H- spin pairs in noe's
    
    for i2=1:9
        fgetl(noe_input_fileID);
    end
    
    while ~feof(noe_input_fileID)
        noe_file=fgetl(noe_input_fileID);
        if noe_file~=-2
            test=strsplit(noe_file,',');
            
            if str2double(test{1,5})==residues(i)
                if strcmp(test{1,6},'H')
                    if str2double(test{1,10})<=max_distance
                        
                        %% atom in spin pair from residue ipeak_list{total,6*(i-1)+3}
                        
                        spin_pair_atom=test{1,9};
                        
                        
                        spin_pair_residue_type=test{1,7};
                        spin_pair_residue=str2double(test{1,8});
                        
                        total=total+1;
                        
                        peak_list{total,6*(i-1)+1}=spin_pair_residue;
                        peak_list{total,6*(i-1)+2}=test{1,9};
                        peak_list{total,6*(i-1)+4}=spin_pair_residue_type;
                        peak_list{total,6*(i-1)+5}=str2double(test{1,12});
                        
                        peak_list{total,6*(i-1)+3}=999;
                        
                        %% find chemical shift
                        
                        frewind(chemical_shiftx2_fileID);
                        
                        if strcmp(spin_pair_atom,'H')==1|strcmp(spin_pair_atom,'HA')==1
                            while ~feof(chemical_shiftx2_fileID)
                                test=fgetl(chemical_shiftx2_fileID);
                                test=fgetl(chemical_shiftx2_fileID);
                                chem_file=fgetl(chemical_shiftx2_fileID);
                                if chem_file~=-1
                                    test=strsplit(chem_file,' ');
                                    
                                    if spin_pair_residue>99
                                        if size(test,2)==8
                                            if spin_pair_residue==str2double(test{1,1})
                                                if strcmp(spin_pair_atom,'H')==1
                                                    peak_list{total,6*(i-1)+3}=str2double(test{1,7});
                                                end
                                                if strcmp(spin_pair_atom,'HA')==1
                                                    peak_list{total,6*(i-1)+3}=str2double(test{1,8});
                                                end
                                                
                                                if strcmp(spin_pair_atom,'HA2')==1
                                                    
                                                    peak_list{total,6*(i-1)+3}=str2double(test{1,8});
                                                    
                                                end
                                                
                                                
                                                if strcmp(spin_pair_atom,'HA3')==1
                                                    
                                                    peak_list{total,6*(i-1)+3}=str2double(test{1,8});
                                                    
                                                end
                                                
                                                if strcmp(spin_pair_atom,'HH')==1
                                                    
                                                    peak_list{total,6*(i-1)+3}=str2double(test{1,7});
                                                    
                                                end
                                                
                                            end
                                        end
                                    end
                                    
                                    if spin_pair_residue<=99
                                        if size(test,2)==9
                                            if spin_pair_residue==str2double(test{1,2})
                                                
                                                if strcmp(spin_pair_atom,'H')==1
                                                    peak_list{total,6*(i-1)+3}=str2double(test{1,8});
                                                end
                                                if strcmp(spin_pair_atom,'HA')==1
                                                    peak_list{total,6*(i-1)+3}=str2double(test{1,9});
                                                end
                                                
                                                if strcmp(spin_pair_atom,'HA2')==1
                                                    
                                                    peak_list{total,6*(i-1)+3}=str2double(test{1,9});
                                                    
                                                end
                                                
                                                
                                                if strcmp(spin_pair_atom,'HA3')==1
                                                    
                                                    peak_list{total,6*(i-1)+3}=str2double(test{1,9});
                                                    
                                                end
                                                
                                                if strcmp(spin_pair_atom,'HH')==1
                                                    
                                                    peak_list{total,6*(i-1)+3}=str2double(test{1,8});
                                                    
                                                end
                                                
                                                
                                            end
                                        end
                                    end
                                    
                                    
                                end
                            end
                        end
                        
                        
                        
                        if strcmp(spin_pair_atom,'H')==0|strcmp(spin_pair_atom,'HA')==0
                            while ~feof(chemical_shiftx2_fileID)
                                chem_file=fgetl(chemical_shiftx2_fileID);
                                if chem_file~=-1
                                    test=strsplit(chem_file,' ');
                                    
                                    if spin_pair_residue>99
                                        if size(test,2)>14
                                            if spin_pair_residue==str2double(test{1,1})
                                                if strcmp(spin_pair_atom,'HB')==1
                                                    peak_list{total,6*(i-1)+3}=str2double(test{1,3});
                                                end
                                                if strcmp(spin_pair_atom,'HB2')==1
                                                    peak_list{total,6*(i-1)+3}=str2double(test{1,4});
                                                end
                                                if strcmp(spin_pair_atom,'HB3')==1
                                                    peak_list{total,6*(i-1)+3}=str2double(test{1,5});
                                                end
                                                if strcmp(spin_pair_atom,'HD1')==1
                                                    peak_list{total,6*(i-1)+3}=str2double(test{1,6});
                                                end
                                                if strcmp(spin_pair_atom,'HD2')==1
                                                    peak_list{total,6*(i-1)+3}=str2double(test{1,7});
                                                end
                                                if strcmp(spin_pair_atom,'HD21')==1
                                                    peak_list{total,6*(i-1)+3}=str2double(test{1,8});
                                                end
                                                if strcmp(spin_pair_atom,'HD22')==1
                                                    peak_list{total,6*(i-1)+3}=str2double(test{1,9});
                                                end
                                                if strcmp(spin_pair_atom,'HD3')==1
                                                    peak_list{total,6*(i-1)+3}=str2double(test{1,10});
                                                end
                                                if strcmp(spin_pair_atom,'HE')==1
                                                    peak_list{total,6*(i-1)+3}=str2double(test{1,11});
                                                end
                                                if strcmp(spin_pair_atom,'HE1')==1
                                                    peak_list{total,6*(i-1)+3}=str2double(test{1,12});
                                                end
                                                if strcmp(spin_pair_atom,'HE2')==1
                                                    peak_list{total,6*(i-1)+3}=str2double(test{1,13});
                                                end
                                                if strcmp(spin_pair_atom,'HE21')==1
                                                    peak_list{total,6*(i-1)+3}=str2double(test{1,14});
                                                end
                                                if strcmp(spin_pair_atom,'HE22')==1
                                                    peak_list{total,6*(i-1)+3}=str2double(test{1,15});
                                                end
                                                if strcmp(spin_pair_atom,'HE3')==1
                                                    peak_list{total,6*(i-1)+3}=str2double(test{1,16});
                                                end
                                                if strcmp(spin_pair_atom,'HG')==1
                                                    peak_list{total,6*(i-1)+3}=str2double(test{1,17});
                                                end
                                                if strcmp(spin_pair_atom,'HG1')==1
                                                    peak_list{total,6*(i-1)+3}=str2double(test{1,18});
                                                end
                                                if strcmp(spin_pair_atom,'HG12')==1
                                                    peak_list{total,6*(i-1)+3}=str2double(test{1,19});
                                                end
                                                if strcmp(spin_pair_atom,'HG13')==1
                                                    peak_list{total,6*(i-1)+3}=str2double(test{1,20});
                                                end
                                                if strcmp(spin_pair_atom,'HG2')==1
                                                    peak_list{total,6*(i-1)+3}=str2double(test{1,21});
                                                end
                                                if strcmp(spin_pair_atom,'HG3')==1
                                                    peak_list{total,6*(i-1)+3}=str2double(test{1,22});
                                                end
                                                if strcmp(spin_pair_atom,'HZ')==1
                                                    peak_list{total,6*(i-1)+3}=str2double(test{1,23});
                                                end
                                                
                                                if strcmp(spin_pair_atom,'HD11')==1|strcmp(spin_pair_atom,'HD12')|strcmp(spin_pair_atom,'HD13')==1
                                                    
                                                    if strcmp(spin_pair_residue_type,"ALA")==1|strcmp(spin_pair_residue_type,"THR")==1 ...
                                                            |strcmp(spin_pair_residue_type,"LEU")==1|strcmp(spin_pair_residue_type,"ILE")==1 ...
                                                            |strcmp(spin_pair_residue_type,"MET")==1|strcmp(spin_pair_residue_type,"VAL")==1
                                                        
                                                        peak_list{total,6*(i-1)+3}=str2double(test{1,6});
                                                        
                                                    end
                                                end
                                                
                                                if strcmp(spin_pair_atom,'HD21')==1|strcmp(spin_pair_atom,'HD22')|strcmp(spin_pair_atom,'HD23')==1
                                                    
                                                    if strcmp(spin_pair_residue_type,"ALA")==1|strcmp(spin_pair_residue_type,"THR")==1 ...
                                                            |strcmp(spin_pair_residue_type,"LEU")==1|strcmp(spin_pair_residue_type,"ILE")==1 ...
                                                            |strcmp(spin_pair_residue_type,"MET")==1|strcmp(spin_pair_residue_type,"VAL")==1
                                                        
                                                        peak_list{total,6*(i-1)+3}=str2double(test{1,7});
                                                        
                                                    end
                                                end
                                                
                                                if strcmp(spin_pair_atom,'HG11')==1
                                                    
                                                    if strcmp(spin_pair_residue_type,"ALA")==1|strcmp(spin_pair_residue_type,"THR")==1 ...
                                                            |strcmp(spin_pair_residue_type,"LEU")==1|strcmp(spin_pair_residue_type,"ILE")==1 ...
                                                            |strcmp(spin_pair_residue_type,"MET")==1|strcmp(spin_pair_residue_type,"VAL")==1
                                                        
                                                        peak_list{total,6*(i-1)+3}=str2double(test{1,18});
                                                        
                                                    end
                                                end
                                                
                                                if strcmp(spin_pair_atom,'HG21')==1|strcmp(spin_pair_atom,'HG22')|strcmp(spin_pair_atom,'HG23')==1
                                                    
                                                    if strcmp(spin_pair_residue_type,"ALA")==1|strcmp(spin_pair_residue_type,"THR")==1 ...
                                                            |strcmp(spin_pair_residue_type,"LEU")==1|strcmp(spin_pair_residue_type,"ILE")==1 ...
                                                            |strcmp(spin_pair_residue_type,"MET")==1|strcmp(spin_pair_residue_type,"VAL")==1
                                                        
                                                        peak_list{total,6*(i-1)+3}=str2double(test{1,21});
                                                        
                                                    end
                                                end
                                                
                                                
                                                if strcmp(spin_pair_atom,'HB1')==1|strcmp(spin_pair_atom,'HB2')|strcmp(spin_pair_atom,'HB3')==1
                                                    
                                                    if strcmp(spin_pair_residue_type,"ALA")==1|strcmp(spin_pair_residue_type,"THR")==1 ...
                                                            |strcmp(spin_pair_residue_type,"LEU")==1|strcmp(spin_pair_residue_type,"ILE")==1 ...
                                                            |strcmp(spin_pair_residue_type,"MET")==1|strcmp(spin_pair_residue_type,"VAL")==1
                                                        
                                                        peak_list{total,6*(i-1)+3}=str2double(test{1,3});
                                                        
                                                    end
                                                end
                                                
                                                
                                            end
                                        end
                                    end
                                    
                                    if spin_pair_residue<=99
                                        if size(test,2)>14
                                            if spin_pair_residue==str2double(test{1,2})
                                                if strcmp(spin_pair_atom,'HB')==1
                                                    peak_list{total,6*(i-1)+3}=str2double(test{1,4});
                                                end
                                                if strcmp(spin_pair_atom,'HB2')==1
                                                    peak_list{total,6*(i-1)+3}=str2double(test{1,5});
                                                end
                                                if strcmp(spin_pair_atom,'HB3')==1
                                                    peak_list{total,6*(i-1)+3}=str2double(test{1,6});
                                                end
                                                if strcmp(spin_pair_atom,'HD1')==1
                                                    peak_list{total,6*(i-1)+3}=str2double(test{1,7});
                                                end
                                                if strcmp(spin_pair_atom,'HD2')==1
                                                    peak_list{total,6*(i-1)+3}=str2double(test{1,8});
                                                end
                                                if strcmp(spin_pair_atom,'HD21')==1
                                                    peak_list{total,6*(i-1)+3}=str2double(test{1,9});
                                                end
                                                if strcmp(spin_pair_atom,'HD22')==1
                                                    peak_list{total,6*(i-1)+3}=str2double(test{1,10});
                                                end
                                                if strcmp(spin_pair_atom,'HD3')==1
                                                    peak_list{total,6*(i-1)+3}=str2double(test{1,11});
                                                end
                                                if strcmp(spin_pair_atom,'HE')==1
                                                    peak_list{total,6*(i-1)+3}=str2double(test{1,12});
                                                end
                                                if strcmp(spin_pair_atom,'HE1')==1
                                                    peak_list{total,6*(i-1)+3}=str2double(test{1,13});
                                                end
                                                if strcmp(spin_pair_atom,'HE2')==1
                                                    peak_list{total,6*(i-1)+3}=str2double(test{1,14});
                                                end
                                                if strcmp(spin_pair_atom,'HE21')==1
                                                    peak_list{total,6*(i-1)+3}=str2double(test{1,15});
                                                end
                                                if strcmp(spin_pair_atom,'HE22')==1
                                                    peak_list{total,6*(i-1)+3}=str2double(test{1,16});
                                                end
                                                if strcmp(spin_pair_atom,'HE3')==1
                                                    peak_list{total,6*(i-1)+3}=str2double(test{1,17});
                                                end
                                                if strcmp(spin_pair_atom,'HG')==1
                                                    peak_list{total,6*(i-1)+3}=str2double(test{1,18});
                                                end
                                                if strcmp(spin_pair_atom,'HG1')==1
                                                    peak_list{total,6*(i-1)+3}=str2double(test{1,19});
                                                end
                                                if strcmp(spin_pair_atom,'HG12')==1
                                                    peak_list{total,6*(i-1)+3}=str2double(test{1,20});
                                                end
                                                if strcmp(spin_pair_atom,'HG13')==1
                                                    peak_list{total,6*(i-1)+3}=str2double(test{1,21});
                                                end
                                                if strcmp(spin_pair_atom,'HG2')==1
                                                    peak_list{total,6*(i-1)+3}=str2double(test{1,22});
                                                end
                                                if strcmp(spin_pair_atom,'HG3')==1
                                                    peak_list{total,6*(i-1)+3}=str2double(test{1,23});
                                                end
                                                if strcmp(spin_pair_atom,'HZ')==1
                                                    peak_list{total,6*(i-1)+3}=str2double(test{1,24});
                                                end
                                            end
                                            
                                            
                                            if strcmp(spin_pair_atom,'HD11')==1|strcmp(spin_pair_atom,'HD12')|strcmp(spin_pair_atom,'HD13')==1
                                                
                                                if strcmp(spin_pair_residue_type,"ALA")==1|strcmp(spin_pair_residue_type,"THR")==1 ...
                                                        |strcmp(spin_pair_residue_type,"LEU")==1|strcmp(spin_pair_residue_type,"ILE")==1 ...
                                                        |strcmp(spin_pair_residue_type,"MET")==1|strcmp(spin_pair_residue_type,"VAL")==1
                                                    
                                                    peak_list{total,6*(i-1)+3}=str2double(test{1,7});
                                                    
                                                end
                                            end
                                            
                                            if strcmp(spin_pair_atom,'HD21')==1|strcmp(spin_pair_atom,'HD22')|strcmp(spin_pair_atom,'HD23')==1
                                                
                                                if strcmp(spin_pair_residue_type,"ALA")==1|strcmp(spin_pair_residue_type,"THR")==1 ...
                                                        |strcmp(spin_pair_residue_type,"LEU")==1|strcmp(spin_pair_residue_type,"ILE")==1 ...
                                                        |strcmp(spin_pair_residue_type,"MET")==1|strcmp(spin_pair_residue_type,"VAL")==1
                                                    
                                                    peak_list{total,6*(i-1)+3}=str2double(test{1,8});
                                                    
                                                end
                                            end
                                            
                                            if strcmp(spin_pair_atom,'HG11')==1
                                                
                                                if strcmp(spin_pair_residue_type,"ALA")==1|strcmp(spin_pair_residue_type,"THR")==1 ...
                                                        |strcmp(spin_pair_residue_type,"LEU")==1|strcmp(spin_pair_residue_type,"ILE")==1 ...
                                                        |strcmp(spin_pair_residue_type,"MET")==1|strcmp(spin_pair_residue_type,"VAL")==1
                                                    
                                                    peak_list{total,6*(i-1)+3}=str2double(test{1,19});
                                                    
                                                end
                                            end
                                            
                                            if strcmp(spin_pair_atom,'HG21')==1|strcmp(spin_pair_atom,'HG22')|strcmp(spin_pair_atom,'HG23')==1
                                                
                                                if strcmp(spin_pair_residue_type,"ALA")==1|strcmp(spin_pair_residue_type,"THR")==1 ...
                                                        |strcmp(spin_pair_residue_type,"LEU")==1|strcmp(spin_pair_residue_type,"ILE")==1 ...
                                                        |strcmp(spin_pair_residue_type,"MET")==1|strcmp(spin_pair_residue_type,"VAL")==1
                                                    
                                                    peak_list{total,6*(i-1)+3}=str2double(test{1,22});
                                                    
                                                end
                                            end
                                            
                                            
                                            if strcmp(spin_pair_atom,'HB1')==1|strcmp(spin_pair_atom,'HB2')|strcmp(spin_pair_atom,'HB3')==1
                                                
                                                if strcmp(spin_pair_residue_type,"ALA")==1|strcmp(spin_pair_residue_type,"THR")==1 ...
                                                        |strcmp(spin_pair_residue_type,"LEU")==1|strcmp(spin_pair_residue_type,"ILE")==1 ...
                                                        |strcmp(spin_pair_residue_type,"MET")==1|strcmp(spin_pair_residue_type,"VAL")==1
                                                    
                                                    peak_list{total,6*(i-1)+3}=str2double(test{1,4});
                                                    
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
        end
    end
end


outputID=fopen(output_file_spectrum,'w');

for i=1:size(peak_list,1)
    for i2=1:size(peak_list,2)/6
        
        fprintf(outputID,'%f,%s,%f,%s,%f,%s,',peak_list{i,(i2-1)*6+1},peak_list{i,(i2-1)*6+2},peak_list{i,(i2-1)*6+3},peak_list{i,(i2-1)*6+4},peak_list{i,(i2-1)*6+5},peak_list{i,(i2-1)*6+6});

    end
    fprintf(outputID,'\n');
end


fclose(noe_input_fileID);
fclose(chemical_shiftx2_fileID);




