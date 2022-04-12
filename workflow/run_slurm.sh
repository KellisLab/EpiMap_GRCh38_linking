#!/usr/bin/env bash
## Derived from https://github.com/rusalkaguy/snakemake-slurm-tutorial


# https://snakemake.readthedocs.io/en/stable/project_info/faq.html#how-do-i-remove-all-files-created-by-snakemake-i-e-like-make-clean
if [ "$1" == "clean" ]; then
    echo 'rm $(snakemake --summary | tail -n+2 | cut -f1)'
    snakemake --summary | tail -n+2 | cut -f1
    rm -f $(snakemake --summary | tail -n+2 | cut -f1)
    exit 0
fi

mkdir -p ../log
SM_PARAMS="job-name ntasks partition time error output"
SM_ARGS="--cpus-per-task 20"
for P in ${SM_PARAMS}; do
    SM_ARGS="$SM_ARGS --$P {cluster.$P}"
done
echo "SM_ARGS: ${SM_ARGS}"


snakemake \
    $* \
     --latency-wait 30 \
    -j 100 \
    --cluster-config ../config/cluster.slurm.luria.json \
    --cluster "sbatch $SM_ARGS" \
    -d .. \
    --use-conda \
    --until all \
    --restart-times 10
