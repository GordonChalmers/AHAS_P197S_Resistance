
%% find ordering of fitnesses 
molecule_fitness=zeros(total_files,1);
for file_idx=1:total_files
    molecule_pdb=readlines("gold_soln_Chlorsulfuron_m1_"+file_idx+".mol");
    molecule_fitness(file_idx)=str2double(molecule_pdb(84));
end
[mol_fit,sort_fit_index]=sort(molecule_fitness,'descend');


total_files=500;
total_atoms=23;
molecule_coords=zeros(total_files,23,3);
for file_idx=1:total_files
    molecule_pdb=readlines("gold_soln_Chlorsulfuron_m1_"+file_idx+".mol");
    for idx=5:4+total_atoms
        text_split=split(molecule_pdb(idx));
        molecule_coords(file_idx,idx-4,:)=[str2double(text_split(2)); str2double(text_split(3)); str2double(text_split(4))];
    end
end

rmsd_molecule=zeros(total_files,1);
centroid=zeros(23,3);
center_file=353;  %% cntr=#7 score and in cavity  P197 
%% center_file=73;  %% S197 cntr=10
diff=zeros(1,3);
for file_idx=1:500
    for idx=1:total_atoms
        diff(:)=molecule_coords(file_idx,idx,:)-molecule_coords(center_file,idx,:);
        rmsd_molecule(file_idx)=rmsd_molecule(file_idx)+diff*diff';
        for q=1:3
            centroid(idx,q)=centroid(idx,q)+molecule_coords(file_idx,idx,q);
        end
    end
end
rmsd_molecule=sqrt(rmsd_molecule/total_atoms);  %% a vector 
centroid=centroid/total_files;

%% average RMSD in order of high to low score, order is sort_fit_index
%% running total from 1 to 500 
cntr=7;
total_rmsd_molecule=zeros(total_files,1);
for idx=cntr:500 
    for idx2=cntr:idx
        total_rmsd_molecule(idx)=total_rmsd_molecule(idx)+rmsd_molecule(sort_fit_index(idx2)); 
    end
    total_rmsd_molecule(idx)=total_rmsd_molecule(idx)/(idx-cntr+1);
end

plot(total_rmsd_molecule(cntr:500))




%% with reference to P197 center rank#7 353
rmsd_molecule=zeros(total_files,1);
centroid=zeros(23,3);
%% center_file=353;  %% #7 score and in cavity  P197 
%% center_file=285;  %% S197
diff=zeros(1,3);
for file_idx=1:500
    for idx=1:total_atoms
        diff(:)=molecule_coords(file_idx,idx,:)-P197_molecule_coords(353,idx,:);
        rmsd_molecule(file_idx)=rmsd_molecule(file_idx)+diff*diff';
        for q=1:3
            centroid(idx,q)=centroid(idx,q)+molecule_coords(file_idx,idx,q);
        end
    end
end
rmsd_molecule=sqrt(rmsd_molecule/total_atoms);  %% a vector 
centroid=centroid/total_files;

%% average RMSD in order of high to low score, order is sort_fit_index
%% running total from 1 to 500 
total_rmsd_molecule=zeros(total_files,1);
for idx=7:500 
    for idx2=7:idx
        total_rmsd_molecule(idx)=total_rmsd_molecule(idx)+rmsd_molecule(sort_fit_index(idx2)); 
    end
    total_rmsd_molecule(idx)=total_rmsd_molecule(idx)/(idx-6);
end

plot(total_rmsd_molecule(7:500))




