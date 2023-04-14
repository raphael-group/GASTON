# sample_values="151507 151508 151509 151510 151669 151670 151671 151672 151673 151674 151675 151676"
sample_values="151669 151670 151671 151672 151673 151674 151675 151676"
trial_values="0 1 2 3 4 5 6 7 8 9"

fS_hidden_values="20 50 100 200"
fA_hidden_values="10"

# optimizer_values="sgd adam"
optimizer_values="adam"

num_cores=28

scripts=/n/fs/ragr-research/projects/network-mutations/manifold-alignment/NN_scripts_ALL
intermediate=/n/fs/ragr-research/projects/network-mutations/manifold-alignment/NN_intermediate_ALL_v2

mkdir -p $intermediate

# Submit jobs.
for sample in $sample_values
do
    for trial in $trial_values
    do
        for nhS in $fS_hidden_values
        do
            for nhA in $fA_hidden_values
            do
                for optimizer in $optimizer_values
                do

                    # Wait before submitting more jobs.
                    while true ; do
                        num_jobs=`squeue -u $USER | wc -l`
                        if [ "$num_jobs" -lt "751" ]; then
                            break
                        fi
                        sleep 30
                    done

                    echo $sample
                    echo $nhS
                    echo $nhA
                    echo $optimizer
                    echo $trial

                    mkdir -p $intermediate/sample_"$sample"_nhS_"$nhS"_nhA_"$nhA"_optimizer_"$optimizer"_trial_"$trial"

                    # NEXT, DO THE NETMIX SCRIPT
                    sbatch -A raphael \
                        -o $intermediate/sample_"$sample"_nhS_"$nhS"_nhA_"$nhA"_optimizer_"$optimizer"_trial_"$trial"/stdout.txt \
                        -e $intermediate/sample_"$sample"_nhS_"$nhS"_nhA_"$nhA"_optimizer_"$optimizer"_trial_"$trial"/stderr.txt \
                        $scripts/./indiv_sbatch.sh \
                            $sample \
                            $nhS \
                            $nhA \
                            $optimizer \
                            $trial 
                    sleep 0.25
                    echo
                done
            done
        done
    done
done