nhS_list="50 100 200"
trial_values="0 1 2 3 4 5 6 7 8 9"
optimizer_list="adam"
partition_list="lower upper"

num_cores=28

scripts=/n/fs/ragr-research/projects/network-mutations/manifold-alignment/olfactory_glmpca/scripts_NN
intermediate=/n/fs/ragr-research/projects/network-mutations/manifold-alignment/olfactory_glmpca/intermediate_NN_v2

mkdir -p $intermediate

# Submit jobs.
for trial in $trial_values
do
    for optimizer in $optimizer_list
    do
        for partition in $partition_list
        do
            for nhS in $nhS_list
            do

                # Wait before submitting more jobs.
                while true ; do
                    num_jobs=`squeue -u $USER | wc -l`
                    if [ "$num_jobs" -lt "751" ]; then
                        break
                    fi
                    sleep 30
                done

                echo $trial
                echo $optimizer
                echo $partition
                echo $nhS

                mkdir -p $intermediate/nhS_"$nhS"_optimizer_"$optimizer"_partition_"$partition"_trial_"$trial"

                # NEXT, DO THE NETMIX SCRIPT
                sbatch -A raphael \
                    -o $intermediate/nhS_"$nhS"_optimizer_"$optimizer"_partition_"$partition"_trial_"$trial"/stdout.txt \
                    -e $intermediate/nhS_"$nhS"_optimizer_"$optimizer"_partition_"$partition"_trial_"$trial"/stderr.txt \
                    $scripts/./indiv_sbatch.sh \
                        $nhS \
                        $optimizer \
                        $partition \
                        $trial 
                sleep 0.25
                echo
            done
        done
    done
done