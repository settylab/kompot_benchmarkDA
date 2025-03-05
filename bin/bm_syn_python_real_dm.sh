#!/bin/bash

module purge
#module load R/4.3.1-gfbf-2022b
module load ImageMagick/7.1.0-53-GCCcore-12.2.0
module load GSL/2.7-GCCcore-12.2.0
module load cuDNN/8.4.1.50-CUDA-11.7.0
eval "$(conda shell.bash hook)"

# set slurm parameters
time=1-00:00:00
partition=campus-new

script_path="$(readlink -f "$0")"
script_dir="$(dirname "$script_path")"
if [ -z "${root+x}" ]; then
    export root="$(readlink -f "$script_dir/..")"
fi
cd ${root}/python_method
echo "files are in :$root/python_method" 

#data_dir="$root/data"
#out_dir="$root/benchmark_python"

#echo "data_dir is $data_dir" 

## Run real data ##

data_id=$1
analysis_layer=$2
iteration_num=$3
balance_bool=$4
mellon_method=$5
norm_or_not=$6
correct=$7
hyper=$8
# job_number=0
# M2 M3 M4 M5 M6 M7
# 44 45
#$(seq 0.75 0.1 0.85)




if [[ "$data_id" == "covid19-pbmc" ]]
    then
    data_dir=${root}/data/real/$data_id
    data_file=/fh/fast/setty_m/user/ryang/differential_abundance/benchmarkDA/data/real/covid19-pbmc/covid_dm.h5ad
    #RBC B CD14_Monocyte CD8_T CD4_T Platelet NK Granulocyte CD16_Monocyte gd_T pDC DC
    pops=$(for m in PB RBC B CD14_Monocyte CD8_T CD4_T Platelet NK Granulocyte CD16_Monocyte gd_T pDC DC; do echo $m; done)
    #R_methods=$(for m in mellon mellon_dm mellon_hls milo daseq cydar cna meld louvain; do echo $m; done)
    R_methods=$(for m in mellon meld cna mellon_high_ls meld_default; do echo $m; done)
    batch_vec=0
    k=30
    resolution=0.5
    beta=25
    downsample=3
    mem=32g
    pop_col="cell.type.coarse"
elif [ "$data_id" == "aging" ]
    then
    data_dir=${root}/data/real/$data_id
    data_file=/fh/fast/setty_m/user/ryang/differential_abundance/benchmarkDA/data/real/aging/aging_dm.h5ad
    #Ery_P HSC ILC Immature_B_cell LMPP MBE MKP Mature_B_cell Mono_P Monocyte Myelo_P NK Neutrophil Pre-B_cell T_cell Treg cDC pDC
    pops=$(for m in CLP Ery_P HSC ILC Immature_B_cell LMPP MBE MKP Mature_B_cell Mono_P Monocyte Myelo_P NK Neutrophil Pre-B_cell T_cell Treg cDC pDC; do echo $m; done)
    #Ery_P HSC ILC Immature_B_cell LMPP MBE MKP Mature_B_cell Mono_P Monocyte Myelo_P NK Neutrophil Pre-B_cell T_cell Treg cDC pDC
    #R_methods=$(for m in mellon mellon_dm mellon_hls milo daseq cydar cna meld louvain milo_batch cna_batch louvain_batch; do echo $m; done)
    #0.75 1 1.25 1.5
    #$(for m in $(seq 0 0.1 1) $(seq 1 0.25 1.75) 20; do echo $m; done)
    R_methods=$(for m in mellon meld cna mellon_high_ls meld_default; do echo $m; done)
    #R_methods=$(for m in mellon meld cna; do echo $m; done)
    #batch_vec=$(for m in $(seq 0 0.1 1); do echo $m; done)
    batch_vec=$(for m in 0.0; do echo $m; done)
    k=30
    resolution=1
    beta=64  
    # beta used ad.X to fit the meld benchmark model and set KNN = 30, find the beta with smallest mse
    # this beta will be changed after the parameter tuning
    downsample=3
    mem=8g
    pop_col="midres_celltype_benchmarking"
elif [[ "$data_id" == "bcr-xl" ]]
    then
    data_dir=${root}/data/real/$data_id
    data_file=/fh/fast/setty_m/user/ryang/differential_abundance/benchmarkDA/data/real/bcr-xl/bcr_xl_anndata_revised.h5ad
    pops=$(for m in CD4_T-cells NK_cells CD8_T-cells B-cells_IgM+ monocytes surface- B-cells_IgM- DC; do echo $m; done)
    R_methods=$(for m in mellon meld cna; do echo $m; done)
    #R_methods=$(for m in mellon mellon_dm mellon_hls milo daseq cydar cna meld louvain; do echo $m; done)
    batch_vec=0
    k=30
    resolution=0.6
    beta=23
    downsample=10
    mem=32g
    pop_col="cell_type"
elif [[ "$data_id" == "pancreas" ]]
    then
    data_dir=${root}/data/real/$data_id
    data_file=/fh/fast/setty_m/user/ryang/differential_abundance/benchmarkDA/data/real/pancreas/pancreas_dm.h5ad
    pops=$(for m in delta_cell alpha_cell gamma_cell acinar_cell beta_cell ductal_cell epsilon_cell; do echo $m; done)
    R_methods=$(for m in mellon meld cna; do echo $m; done)
    #R_methods=$(for m in milo daseq cydar cna meld louvain; do echo $m; done)
    batch_vec=0
    k=30
    resolution=1.2
    beta=80
    downsample=3
    mem=8g
    pop_col="cell_type"
elif [[ "$data_id" == "levine32" ]]
    then
    data_dir=${root}/data/real/$data_id
    data_file=/fh/fast/setty_m/user/ryang/differential_abundance/benchmarkDA/data/real/levine32/levine32_anndata_revised_deduplicated.h5ad
    # CD4_T_cells CD8_T_cells Pre_B_cells Mature_B_cells Monocytes Basophils
    pops=$(for m in pDCs CD4_T_cells CD8_T_cells Pre_B_cells Mature_B_cells Monocytes Basophils; do echo $m; done)
    R_methods=$(for m in mellon meld cna; do echo $m; done)
    #R_methods=$(for m in milo daseq cydar cna meld louvain; do echo $m; done)
    batch_vec=0
    k=30
    resolution=0.6
    beta=36
    downsample=25
    mem=96g
    pop_col="cell_type"
fi


#PB CD14_Monocyte CD8_T CD4_T Platelet NK Granulocyte CD16_Monocyte gd_T pDC DC
echo "data_dir is $data_dir" 
job_number=0


for p in $pops;
    do
    # 44 45
    for seed in 43 44 45
        do
        for enr in $(seq 0.75 0.1 0.95)
        	do
			for batch_sd_num in $batch_vec
                do
                for method in $R_methods
                    do
                    for iteration in $(seq 0 1 $iteration_num)
                        do
                            ((job_number++))
                            if [ -z "$SLURM_ARRAY_TASK_ID" ] || [ "$job_number" -ne "$SLURM_ARRAY_TASK_ID" ]; then
                                continue
                            fi
                            jobid=${data_id}-${p}-${enr}-${seed}-${batch_sd_num}-${balance_bool}-${analysis_layer}
                            jobid_2=${data_id}-${p}-${enr}-${seed}-${batch_sd_num}-${balance_bool}-${analysis_layer}-${mellon_method}-${norm_or_not}-${hyper}-${correct}

                            save_path=${root}/benchmark_dm/real/$data_id/${jobid_2}
                            save_path_iteration=${root}/benchmark_dm/real/$data_id/${jobid_2}/iteration_${iteration}
                            mkdir -p "$save_path_iteration"
                            echo "Doing $jobid ..."
                            if [[ "$method" == "mellon" ]]; then

                                conda deactivate
                                conda activate DiffAbundance
                                python Mellon_bm.py \
                                    --file_path ${data_file} \
                                    --pop ${p} \
                                    --pop_enr $enr \
                                    --pop_column ${pop_col} \
                                    --ds_type $data_id \
                                    --batch_sd ${batch_sd_num} \
                                    --input_file $data_dir/${jobid}/ \
                                    --package $method \
                                    --seed ${seed} \
                                    --layer_embedding X_pca \
                                    --n_dm 0 \
                                    --mellon_d_method ${mellon_method} \
                                    --norm_density ${norm_or_not} \
                                    --hyperparameter ${hyper} \
                                    --corrected ${corrected} \
                                    --ls_factor 1.5 \
                                    --ls_mode PCA \
                                    --output_dir $save_path/
                                exit $!
                            elif [[ "$method" == "mellon_high_ls" ]]; then
                                conda deactivate
                                conda activate DiffAbundance
                                python Mellon_bm.py \
                                    --file_path ${data_file} \
                                    --pop ${p} \
                                    --pop_enr $enr \
                                    --pop_column ${pop_col} \
                                    --ds_type $data_id \
                                    --batch_sd ${batch_sd_num} \
                                    --input_file $data_dir/${jobid}/ \
                                    --package $method \
                                    --seed ${seed} \
                                    --layer_embedding X_pca \
                                    --n_dm 0 \
                                    --mellon_d_method ${mellon_method} \
                                    --norm_density ${norm_or_not} \
                                    --hyperparameter ${hyper} \
                                    --corrected ${corrected} \
                                    --ls_factor 10 \
                                    --ls_mode DM \
                                    --output_dir $save_path/
                                exit $!
                            elif [[ "$method" == "meld" ]]; then
                                conda deactivate
                                conda activate DiffAbundance
                                python meld_bm.py \
                                    --file_path ${data_file} \
                                    --pop ${p} \
                                    --pop_enr $enr \
                                    --pop_column ${pop_col} \
                                    --ds_type $data_id \
                                    --batch_sd ${batch_sd_num} \
                                    --input_file $data_dir/${jobid}/ \
                                    --package $method \
                                    --seed ${seed} \
                                    --layer_embedding X_pca \
                                    --beta $beta \
                                    --k_meld $k \
                                    --output_dir $save_path/
                                exit $!
                            elif [[ "$method" == "meld_default" ]]; then
                                conda deactivate
                                conda activate DiffAbundance
                                python meld_bm.py \
                                    --file_path ${data_file} \
                                    --pop ${p} \
                                    --pop_enr $enr \
                                    --pop_column ${pop_col} \
                                    --ds_type $data_id \
                                    --batch_sd ${batch_sd_num} \
                                    --input_file $data_dir/${jobid}/ \
                                    --package $method \
                                    --seed ${seed} \
                                    --layer_embedding X_pca \
                                    --beta 40 \
                                    --k_meld $k \
                                    --output_dir $save_path/
                                exit $!
                            elif [[ "$method" == "cna" ]]; then
                                conda deactivate
                                conda activate DiffAbundance
                                python CNA_bm.py \
                                    --file_path ${data_file} \
                                    --pop ${p} \
                                    --pop_enr $enr \
                                    --pop_column ${pop_col} \
                                    --ds_type $data_id \
                                    --batch_sd ${batch_sd_num} \
                                    --input_file $data_dir/${jobid}/ \
                                    --package $method \
                                    --seed ${seed} \
                                    --layer_embedding X_pca \
                                    --k_cna $k \
                                    --output_dir $save_path/
                                exit $!
                            fi                     
                    done					
				done
			done
        done
    done
done

# Submit a slurm array job
jobid="mellon_syn_real_$data_id"
cmd="sbatch -J '$jobid' --time=$time --partition=$partition \
--mem 8g --out '$root/SlurmLog/${jobid_2}_%N_%A_%a.out' --array=1-$job_number \
'$script_path' $1 $2 $3 $4 $5 $6 $7 $8"
echo "$cmd"
eval "$cmd"