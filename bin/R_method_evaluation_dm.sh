#!/bin/bash

# add slurm module first
module purge
module load fhR/4.3.1-foss-2022b
module load ImageMagick/7.1.0-53-GCCcore-12.2.0
module load GSL/2.7-GCCcore-12.2.0
module load cuDNN/8.4.1.50-CUDA-11.7.0
eval "$(conda shell.bash hook)"

# set slurm parameters
time=1-00:00:00
partition=campus-new

data_id=$1
analysis_layer=$2
iteration_num=$3
balance_bool=$4
mellon_method=$5
norm_or_not=$6
correct=$7
hyper=$8

script_path="$(readlink -f "$0")"
script_dir="$(dirname "$script_path")"
if [ -z "${root+x}" ]; then
    export root="$(readlink -f "$script_dir/..")"
fi
cd ${root}



#### NOTE : remember to change the out_dir to make sure that they are in the ideal path

if [ "$data_id" == "cluster" ]
    then
    data_dir=${root}/data/synthetic/$data_id
    data_file=/fh/fast/setty_m/user/ryang/differential_abundance/manuscript_preparation/cluster_dm_10.rds
    pops=$(for p in $(seq 1 1 3); do echo M$p; done)
    R_methods=$(for m in milo daseq cydar louvain; do echo $m; done)
    batch_vec=$(for m in 0 0.75 1 1.25 1.5; do echo $m; done)
    k=30
    resolution=0.2
    beta=33
    downsample=3
    mem=8g
    pop_col="celltype"
    out_dir=${root}/benchmark_dm/synthetic/$data_id
# elif [ "$data_id" == "cluster_balanced" ]
#     then
#     data_dir=${root}/data/synthetic/$data_id
#     data_file=${data_dir}/${data_id}_data_bm.RDS
#     pops=$(for p in $(seq 1 1 3); do echo M$p; done)
#     R_methods=$(for m in milo daseq cydar louvain; do echo $m; done)
#     batch_vec=$(for m in 0 0.75 1 1.25 1.5; do echo $m; done)
#     k=30
#     resolution=0.2
#     beta=33
#     downsample=3
#     mem=8g
#     pop_col="celltype"
#     out_dir=${root}/benchmark_python/synthetic/$data_id
elif [ "$data_id" == "linear" ]
    then
    data_dir=${root}/data/synthetic/$data_id
    data_file=/fh/fast/setty_m/user/ryang/differential_abundance/manuscript_preparation/linear_dm_10.rds
    pops=$(for p in $(seq 1 1 7); do echo M$p; done)
    R_methods=$(for m in milo daseq cydar louvain; do echo $m; done)
    #0.75 1 1.25 1.5
    batch_vec=$(for m in 0 0.75 1 1.25 1.5; do echo $m; done)
    k=30
    resolution=1
    beta=71
    downsample=3
    mem=8g
    pop_col="celltype"
    out_dir=${root}/benchmark_dm/synthetic/$data_id
elif [ "$data_id" == "branch" ]
    then
    data_dir=${root}/data/synthetic/$data_id
    data_file=/fh/fast/setty_m/user/ryang/differential_abundance/manuscript_preparation/branch_dm_10.rds
    pops=$(for p in $(seq 1 1 8); do echo M$p; done)
    R_methods=$(for m in milo daseq cydar louvain; do echo $m; done)
    batch_vec=$(for m in 0 0.75 1 1.25 1.5; do echo $m; done)
    k=30
    resolution=1
    beta=65
    downsample=3
    mem=8g
    pop_col="celltype"
    out_dir=${root}/benchmark_dm/synthetic/$data_id
elif [ "$data_id" == "aging" ]
    then
    data_dir=${root}/data/real/$data_id
    data_file=/fh/fast/setty_m/user/ryang/differential_abundance/manuscript_preparation/aging_dm_30.rds
    pops=$(for m in CLP Ery_P HSC ILC Immature_B_cell LMPP MBE MKP Mature_B_cell Mono_P Monocyte Myelo_P NK Neutrophil Pre-B_cell T_cell Treg cDC pDC; do echo $m; done)
    #Ery_P HSC ILC Immature_B_cell LMPP MBE MKP Mature_B_cell Mono_P Monocyte Myelo_P NK Neutrophil Pre-B_cell T_cell Treg cDC pDC
    #R_methods=$(for m in mellon mellon_dm mellon_hls milo daseq cydar cna meld louvain milo_batch cna_batch louvain_batch; do echo $m; done)
    #0.75 1 1.25 1.5
    #$(for m in $(seq 0 0.1 1) $(seq 1 0.25 1.75) 20; do echo $m; done)
    R_methods=$(for m in milo daseq cydar louvain; do echo $m; done)
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
    out_dir=${root}/benchmark_dm/real/$data_id
elif [[ "$data_id" == "covid19-pbmc" ]]
    then
    data_dir=${root}/data/real/$data_id
    data_file=/fh/fast/setty_m/user/ryang/differential_abundance/manuscript_preparation/covid_dm.rds
    #RBC B
    pops=$(for m in RBC B PB CD14_Monocyte CD8_T CD4_T Platelet NK Granulocyte CD16_Monocyte gd_T pDC DC; do echo $m; done)
    #R_methods=$(for m in mellon mellon_dm mellon_hls milo daseq cydar cna meld louvain; do echo $m; done)
    R_methods=$(for m in milo daseq cydar louvain; do echo $m; done)
    batch_vec=0
    k=30
    resolution=0.5
    beta=25
    downsample=3
    mem=32g
    pop_col="cell.type.coarse"
    out_dir=${root}/benchmark_dm/real/$data_id
elif [[ "$data_id" == "bcr-xl" ]]
    then
    data_dir=${root}/data/real/$data_id
    data_file=/fh/fast/setty_m/user/ryang/differential_abundance/benchmarkDA/data/real/bcr-xl/bcr_xl_preprocessed_sce.rds
    pops=$(for m in CD4_T-cells NK_cells CD8_T-cells B-cells_IgM+ monocytes surface- B-cells_IgM- DC; do echo $m; done)
    R_methods=$(for m in milo daseq cydar louvain; do echo $m; done)
    #R_methods=$(for m in mellon mellon_dm mellon_hls milo daseq cydar cna meld louvain; do echo $m; done)
    batch_vec=0
    k=30
    resolution=0.6
    beta=23
    downsample=10
    mem=32g
    pop_col="cell_type"
    out_dir=${root}/benchmark_dm/real/$data_id
elif [[ "$data_id" == "pancreas" ]]
    then
    data_dir=${root}/data/real/$data_id
    data_file=/fh/fast/setty_m/user/ryang/differential_abundance/manuscript_preparation/pancreas_dm_30.rds
    pops=$(for m in delta_cell alpha_cell gamma_cell acinar_cell beta_cell ductal_cell epsilon_cell; do echo $m; done)
    R_methods=$(for m in milo daseq cydar louvain; do echo $m; done)
    #R_methods=$(for m in milo daseq cydar cna meld louvain; do echo $m; done)
    batch_vec=0
    k=30
    resolution=1.2
    beta=80
    downsample=3
    mem=8g
    pop_col="Factor.Value.inferred.cell.type...authors.labels."
    out_dir=${root}/benchmark_dm/real/$data_id
elif [[ "$data_id" == "levine32" ]]
    then
    data_dir=${root}/data/real/$data_id
    data_file=/fh/fast/setty_m/user/ryang/differential_abundance/benchmarkDA/data/real/levine32/deduplicated_levine32.rds
    # CD4_T_cells CD8_T_cells Pre_B_cells Mature_B_cells Monocytes Basophils
    pops=$(for m in pDCs CD4_T_cells CD8_T_cells Pre_B_cells Mature_B_cells Monocytes Basophils; do echo $m; done)
    R_methods=$(for m in milo daseq cydar louvain; do echo $m; done)
    #R_methods=$(for m in milo daseq cydar cna meld louvain; do echo $m; done)
    batch_vec=0
    k=30
    resolution=0.6
    beta=36
    downsample=25
    mem=96g
    pop_col="cell_type"
    out_dir=${root}/benchmark_dm/real/$data_id
fi


echo "going to start, $data_id"


job_number=0
## Run
for pop in $pops
    do
    echo "the population is $pop"
    for pop_enr in $(seq 0.75 0.1 0.95)
        do
        echo "$pop_enr"
        for seed in 43 44 45
            do
            echo "$seed"
            for batch_sd in $batch_vec
                do
                echo "$batch_sd"
                for iteration in $(seq 0 1 $iteration_num)
                    do
                    for method in $R_methods
                        do
                        echo "$method"
                        ((job_number++))
                        if [ -z "$SLURM_ARRAY_TASK_ID" ] || [ "$job_number" -ne "$SLURM_ARRAY_TASK_ID" ]; then
                            continue
                        fi
                        jobid_2=${data_id}-${pop}-${pop_enr}-${seed}-${batch_sd}-${balance_bool}-${analysis_layer}-${mellon_method}-${norm_or_not}-${hyper}-${correct}

                        ## this jobid_old needs to match with the path where the synthetic labals generated from "modified_benchmarkda_dm_all.sh" are saved 
                        jobid_old=${data_id}-${pop}-${pop_enr}-${seed}-${batch_sd}-${balance_bool}-${analysis_layer}-DM

                        save_path_iteration=${out_dir}/${jobid_2}/iteration_${iteration}/
                        mkdir -p "$save_path_iteration"
                            echo "Doing $jobid ..."
                        echo "starting"
                        Rscript scripts/run_DA_test.r \
                        ## tol parameter for cydar needs to be regulated inside the R script
                            ${data_file} $method $seed $pop \
                            --data_dir ${data_dir}/${jobid_old}/iteration_${iteration}/ \
                            --pop_enrichment $pop_enr \
                            --data_id $data_id \
                            --k $k \
                            --resolution ${resolution} \
                            --downsample ${downsample} \
                            --batchEffect_sd $batch_sd \
                            --outdir ${out_dir}/${jobid_2}/iteration_${iteration}/
                        exit $!
                    done
                done
            done
        done
    done
done

echo $job_number

# Submit a slurm array job
jobid="benchmarkDA_syn_real_$data_id"
cmd="sbatch -J '$jobid' --time=$time --partition=$partition \
--mem $mem --out '$root/SlurmLog/${jobid}_%N_%A_%a.out' --array=1-$job_number \
'$script_path' $1 $2 $3 $4 $5 $6 $7 $8"
echo "$cmd"
eval "$cmd"