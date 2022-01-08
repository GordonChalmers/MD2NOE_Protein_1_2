
%% This program creates an noe spectrum from the 1/r^6 calculation.

%% Inputs are a ppm_one file of chemical shifts, 1/r^6
%%  calculations from the SingleFrameNMR program, and an output
%%  file.

%% The output is a peak list which can be used in noe_vectors_from_peak_list.m.

chemical_shift_ppm_one_fileID=fopen(input_chemical_ppm_one_file,'r');
chemical_shift_ppm_one_second_fileID=fopen(input_chemical_ppm_one_second_file,'r');
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

            if str2double(test{1,3})==residues(i)
                if strcmp(test{1,4},'H')
                    if str2double(test{1,8})<=max_distance

                        %% atom in spin pair from residue i

                        spin_pair_atom=test{1,7};

                        spin_pair_residue_type=test{1,5};
                        spin_pair_residue=str2double(test{1,6});

noe_file=test{1,10};

                        total=total+1;

                        peak_list{total,6*(i-1)+1}=spin_pair_residue;
                        peak_list{total,6*(i-1)+2}=spin_pair_atom;
                        peak_list{total,6*(i-1)+4}=spin_pair_residue_type;
%% noe intensity        peak_list{total,6*(i-1)+5}=str2double(test{1,12});

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


                        frewind(chemical_shift_ppm_one_second_fileID);
                        chem_file=fgetl(chemical_shift_ppm_one_second_fileID);

                        while ~feof(chemical_shift_ppm_one_second_fileID)
                            chem_file=fgetl(chemical_shift_ppm_one_second_fileID);
                            if chem_file~=-1
                                test=strsplit(chem_file,' ');

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


                        noe_fileID=fopen(noe_file,'r');

                       while ~feof(noe_fileID)
                             mixing_noe=fgetl(noe_fileID);
                             if mixing_noe~=-1
                                test=strsplit(mixing_noe,' ');

                                if test{1,1}=mixing_time
                                   noe=str2double(test{1,2});
                                end
                             end
                       end

                       peak_list{total,6*(i-1)+5}=noe;



                    end
                end
            end
        end
    end
end


fclose(noe_input_fileID);
fclose(chemical_shift_ppm_one_fileID);


xlswrite(output_file_spectrum,cell2mat(peak_list));
