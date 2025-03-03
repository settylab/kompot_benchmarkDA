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

#daseq cydar louvain milo_batch cydar_batch louvain_batch


if [ "$data_id" == "cluster" ]
    then
    data_dir=${root}/data/synthetic/$data_id
    data_file=/fh/fast/setty_m/user/ryang/differential_abundance/benchmarkDA/data/synthetic/cluster/cluster_data_bm.RDS
    pops=$(for p in $(seq 1 1 3); do echo M$p; done)
    R_methods=$(for m in milo daseq cydar louvain; do echo $m; done)
    batch_vec=$(for m in 0 0.75 1 1.25 1.5; do echo $m; done)
    k=30
    resolution=0.2
    beta=33
    downsample=3
    mem=8g
    pop_col="celltype"
    out_dir=${root}/benchmark_python/synthetic/$data_id
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
    data_file=/fh/fast/setty_m/user/ryang/differential_abundance/benchmarkDA/data/synthetic/linear/linear_data_bm.RDS
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
    out_dir=${root}/benchmark_python/synthetic/$data_id
elif [ "$data_id" == "branch" ]
    then
    data_dir=${root}/data/synthetic/$data_id
    data_file=/fh/fast/setty_m/user/ryang/differential_abundance/benchmarkDA/data/synthetic/branch/branch_data_bm.RDS
    pops=$(for p in $(seq 1 1 8); do echo M$p; done)
    R_methods=$(for m in milo daseq cydar louvain; do echo $m; done)
    batch_vec=$(for m in 0 0.75 1 1.25 1.5; do echo $m; done)
    k=30
    resolution=1
    beta=65
    downsample=3
    mem=8g
    pop_col="celltype"
    out_dir=${root}/benchmark_python/synthetic/$data_id
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

                        jobid_old=${data_id}-${pop}-${pop_enr}-${seed}-${batch_sd}-${balance_bool}-${analysis_layer}
                        echo "Doing $jobid ..."

                        save_path_iteration=${out_dir}/${jobid_2}/iteration_${iteration}/
                        mkdir -p "$save_path_iteration"
                            echo "Doing $jobid ..."
                        echo "starting"
                        
                        Rscript scripts/run_DA_test.r \
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
--mem $mem --out '$root/SlurmLog_R_test_0624/${jobid}_%N_%A_%a.out' --array=1-$job_number \
'$script_path' $1 $2 $3 $4 $5 $6 $7 $8"
echo "$cmd"
eval "$cmd"