

%% The output is a peak list which can be used in noe_vectors_from_peak_list.m.

chemical_shiftx2_fileID=fopen(input_chemical_shiftx2_file,'r');

mixing=mixing_milliseconds;
residues=input_residues;
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
            
            if str2double(test{1,1})==residues(i)
                
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
                    
                endhttps://outlook.office.com/owa/?realm=outlook.uga.edu&path=/mail/inbox
                
                fclose(noe_fileID);
                
                peak_list{total,6*(i-1)+5}=noe;
                
                peak_list{total,6*(i-1)+3}=999;
                
                
                
                %% Inputs are a shiftx2 file of chemical shifts
                
                %% The output is a peak list which can be used in noe_vectors_from_peak_lhttps://outlook.office.com/owa/?realm=outlook.uga.edu&path=/mail/inboxist.m.
                
                
                frewind(chemical_shiftx2_fileID);
                
                if strcmp(spin_pair_atom,'H')==1|strcmp(spin_pair_atom,'HA')==1
                    while ~feof(chemical_shiftx2_fileID)
                        chem_file=fgetl(chemical_shiftx2_fileID);
                        if chem_file~=-1
                            test=strsplit(chem_file,' ');
                            
                            
                            if spin_pair_residue==str2double(test{1,1})
                                
                                if strcmp(spin_pair_atom,'HA')==1
                                    peak_list{total,6*(i-1)+3}=str2double(test{1,5});
                                end
                                
                                if strcmp(spin_pair_atom,'HA2')==1
                                    
                                    peak_list{total,6*(i-1)+3}=str2double(test{1,5});
                                    
                                end
                                
                                
                                if strcmp(spin_pair_atom,'HA3')==1
                                    
                                    peak_list{total,6*(i-1)+3}=str2double(test{1,5});
                                    
                                end
                                
                                if strcmp(spin_pair_atom,'HH')==1
                                    
                                    peak_list{total,6*(i-1)+3}=str2double(test{1,5});
                                    
                                end
                                
                            end
                            
                            
                        end
                    end
                    
                end
                
                
                
                
                while ~feof(chemical_shiftx2_fileID)
                    
                    chem_file=fgetl(chemical_shiftx2_fileID);
                    if chem_file~=-1
                        test=strsplit(chem_file,' ');
                        
                        if spin_pair_residue==str2double(test{1,1})
                            if strcmp(spin_pair_atom,'HB')==1
                                peak_list{total,6*(i-1)+3}=str2double(test{1,5});
                            end
                            if strcmp(spin_pair_atom,'HB2')==1
                                peak_list{total,6*(i-1)+3}=str2double(test{1,5});
                            end
                            if strcmp(spin_pair_atom,'HB3')==1
                                peak_list{total,6*(i-1)+3}=str2double(test{1,5});
                            end
                            if strcmp(spin_pair_atom,'HD1')==1
                                peak_list{total,6*(i-1)+3}=str2double(test{1,5});
                            end
                            if strcmp(spin_pair_atom,'HD2')==1
                                peak_list{total,6*(i-1)+3}=str2double(test{1,5});
                            end
                            if strcmp(spin_pair_atom,'HD21')==1
                                peak_list{total,6*(i-1)+3}=str2double(test{1,5});
                            end
                            if strcmp(spin_pair_atom,'HD22')==1
                                peak_list{total,6*(i-1)+3}=str2double(test{1,5});
                            end
                            if strcmp(spin_pair_atom,'HD3')==1
                                peak_list{total,6*(i-1)+3}=str2double(test{1,5});
                            end
                            if strcmp(spin_pair_atom,'HE')==1
                                peak_list{total,6*(i-1)+3}=str2double(test{1,5});
                            end
                            if strcmp(spin_pair_atom,'HE1')==1
                                peak_list{total,6*(i-1)+3}=str2double(test{1,5});
                            end
                            if strcmp(spin_pair_atom,'HE2')==1
                                peak_list{total,6*(i-1)+3}=str2double(test{1,5});
                            end
                            if strcmp(spin_pair_atom,'HE21')==1
                                peak_list{total,6*(i-1)+3}=str2double(test{1,5});
                            end
                            if strcmp(spin_pair_atom,'HE22')==1
                                peak_list{total,6*(i-1)+3}=str2double(test{1,5});
                            end
                            if strcmp(spin_pair_atom,'HE3')==1
                                peak_list{total,6*(i-1)+3}=str2double(test{1,5});
                            end
                            if strcmp(spin_pair_atom,'HG')==1
                                peak_list{total,6*(i-1)+3}=str2double(test{1,5});
                            end
                            if strcmp(spin_pair_atom,'HG1')==1
                                peak_list{total,6*(i-1)+3}=str2double(test{1,5});
                            end
                            if strcmp(spin_pair_atom,'HG12')==1
                                peak_list{total,6*(i-1)+3}=str2double(test{1,5});
                            end
                            if strcmp(spin_pair_atom,'HG13')==1
                                peak_list{total,6*(i-1)+3}=str2double(test{1,5});
                            end
                            if strcmp(spin_pair_atom,'HG2')==1
                                peak_list{total,6*(i-1)+3}=str2double(test{1,5});
                            end
                            if strcmp(spin_pair_atom,'HG3')==1
                                peak_list{total,6*(i-1)+3}=str2double(test{1,5});
                            end
                            if strcmp(spin_pair_atom,'HZ')==1
                                peak_list{total,6*(i-1)+3}=str2double(test{1,5});
                            end
                            
                            if strcmp(spin_pair_atom,'HD11')==1|strcmp(spin_pair_atom,'HD12')|strcmp(spin_pair_atom,'HD13')==1
                                
                                if strcmp(spin_pair_residue_type,"ALA")==1|strcmp(spin_pair_residue_type,"THR")==1 ...
                                        |strcmp(spin_pair_residue_type,"LEU")==1|strcmp(spin_pair_residue_type,"ILE")==1 ...
                                        |strcmp(spin_pair_residue_type,"MET")==1|strcmp(spin_pair_residue_type,"VAL")==1
                                    
                                    peak_list{total,6*(i-1)+3}=str2double(test{1,5});
                                    
                                end
                            end
                            
                            if strcmp(spin_pair_atom,'HD21')==1|strcmp(spin_pair_atom,'HD22')|strcmp(spin_pair_atom,'HD23')==1
                                
                                if strcmp(spin_pair_residue_type,"ALA")==1|strcmp(spin_pair_residue_type,"THR")==1 ...
                                        |strcmp(spin_pair_residue_type,"LEU")==1|strcmp(spin_pair_residue_type,"ILE")==1 ...
                                        |strcmp(spin_pair_residue_type,"MET")==1|strcmp(spin_pair_residue_type,"VAL")==1
                                    
                                    peak_list{total,6*(i-1)+3}=str2double(test{1,5});
                                    
                                end
                            end
                            
                            if strcmp(spin_pair_atom,'HG11')==1
                                
                                if strcmp(spin_pair_residue_type,"ALA")==1|strcmp(spin_pair_residue_type,"THR")==1 ...
                                        |strcmp(spin_pair_residue_type,"LEU")==1|strcmp(spin_pair_residue_type,"ILE")==1 ...
                                        |strcmp(spin_pair_residue_type,"MET")==1|strcmp(spin_pair_residue_type,"VAL")==1
                                    
                                    peak_list{total,6*(i-1)+3}=str2double(test{1,5});
                                    
                                end
                            end
                            
                            if strcmp(spin_pair_atom,'HG21')==1|strcmp(spin_pair_atom,'HG22')|strcmp(spin_pair_atom,'HG23')==1
                                
                                if strcmp(spin_pair_residue_type,"ALA")==1|strcmp(spin_pair_residue_type,"THR")==1 ...
                                        |strcmp(spin_pair_residue_type,"LEU")==1|strcmp(spin_pair_residue_type,"ILE")==1 ...
                                        |strcmp(spin_pair_residue_type,"MET")==1|strcmp(spin_pair_residue_type,"VAL")==1
                                    
                                    peak_list{total,6*(i-1)+3}=str2double(test{1,51});
                                    
                                end
                            end
                            
                            
                            if strcmp(spin_pair_atom,'HB1')==1|strcmp(spin_pair_atom,'HB2')|strcmp(spin_pair_atom,'HB3')==1
                                
                                if strcmp(spin_pair_residue_type,"ALA")==1|strcmp(spin_pair_residue_type,"THR")==1 ...
                                        |strcmp(spin_pair_residue_type,"LEU")==1|strcmp(spin_pair_residue_type,"ILE")==1 ...
                                        |strcmp(spin_pair_residue_type,"MET")==1|strcmp(spin_pair_residue_type,"VAL")==1
                                    
                                    peak_list{total,6*(i-1)+3}=str2double(test{1,5});
                                    
                                end
                            end
                            
                            
                        end
                    end
                end
                
                
                
            end
        end
    end
end




fclose(noe_input_fileID);
fclose(chemical_shiftx2_fileID);






