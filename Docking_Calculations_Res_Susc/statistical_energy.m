
%%%%%%%%%%%%%   
%%
%%  Calculate free energy from scores 

free_energy=zeros(total_molecules,total_stereoisomers);
int_energy=zeros(total_molecules,total_stereoisomers);
entropic_energy=zeros(total_molecules,total_stereoisomers);

for i=1:total_molecules 
    for j=1:total_stereoisomers
        
molecule_number_in_bins=number_in_bins(i,j,:);        
molecule_bin_values=bin_values(i,j,:);

T=298;
beta_constant=1.86462*(298/T);  %% 1/(k_B T)
gamma_constant=1/6.5;

%% N\left(\beta\right)=\sum_{i}{N(Z_i)e^{-\beta E(Z_i)}}
%% normalization
N_beta=0;
for bin=1:num_bins 
    N_beta=N_beta+molecule_number_in_bins(bin)*exp(beta_constant*molecule_bin_values(bin)*gamma_constant);
end
%% check at T->large

%% state probabilities
p_probability=zeros(num_bins,1);
for bin=1:num_bins
    p_probability(bin)=molecule_number_in_bins(bin)/N_beta*exp(beta_constant*molecule_bin_values(bin)*gamma_constant);
end

%% interaction energy 
exp_energy=0;
for bin=1:num_bins
    exp_energy=exp_energy+p_probability(bin)*molecule_bin_values(bin)*gamma_constant;
end
int_energy(i,j)=exp_energy;

%% entropic energy 
s_energy=0;
for bin=1:num_bins 
    if p_probability(bin)~=0
        s_energy=s_energy-(p_probability(bin)*log(p_probability(bin)))/beta_constant;
    end
end
entropic_energy(i,j)=s_energy;
    
free_energy(i,j)=int_energy(i,j)+entropic_energy(i,j);

    end
end
