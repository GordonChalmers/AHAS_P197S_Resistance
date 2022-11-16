

total_files=500;
center_file=108;
molecule_coords=zeros(3,total_files);
for file_idx=1:total_files
    molecule_pdb=readlines("gold_soln_Chlorsulfuron_m1_"+file_idx+".mol");
    for idx=5:22
        text_split=split(molecule_pdb(idx));
        molecule_coords(file_idx,idx,:)=[str2double(text_split(2)) str2double(text_split(3)) str2double(text_split(4))];
    end
end

rmsd=0;
centroid=zeros(3,1);
for file_idx=1:500
    rmsd=rmsd+norm(molecule_coords(file_idx,idx,:)-molecule_coords(center_file,idx,:));
    centroid=centroid+molecule_coords(file_idx,:);
end

rmsd=rmsd/(total-1);
centroid=centroid/total;



