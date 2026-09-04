#!/bin/bash -f

if [ $# -lt 1 ]; then
    echo "Usage: $0 <case_name>"
    exit 1
fi

case_name="$1"  # First argument is the case name
num_params="${2:-50}"  # Second argument is the number of batches to do, default 50 
N_it="${3:-1000}" # Third argument is the number of runs to do per batch, default 1000

python_script="../run_Venus_evolution.py"
concat_script="./concatenate_sorting.py"
INPUT_FILE="./venus_run_input"
# Read one argument per line
mapfile -t RAW_ARGS < $INPUT_FILE


# Process and evaluate expressions safely (integers and floats)
ARGS=()
ARGS+=("$N_it")
for arg in "${RAW_ARGS[@]}"; do
    if [[ "$arg" =~ [\*\+/] ]]; then
        # Evaluate with bc (handles floats too)
        result=$(echo "$arg" | bc -l)
        ARGS+=("$result")
    else
        ARGS+=("$arg")
    fi
done


echo "Running with args: ${ARGS[*]}"
# Array to store job IDs
job_ids=()

# Loop through each parameter value from 1 to num_params
for param in $(seq 1 $num_params); do
    dir_name="num_$param"
    mkdir -p "$dir_name"
    cd "$dir_name"

    job_name="num_${param}_${case_name}.sh"

    python_cmd="python3 $python_script"
    for arg in "${ARGS[@]}"; do
        python_cmd+=" \"$arg\""
    done

    # Create SLURM job script for the parameter
    cat <<EOF > "$job_name"   
#!/bin/bash
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -t 60
#SBATCH -p sapphire
#SBATCH --mem-per-cpu=5000
#SBATCH -o run_Venus_evolve.out
#SBATCH -e run_Venus_evolve.err

$python_cmd
EOF

    # Submit the job and store its job ID
    job_id=$(sbatch "$job_name" | awk '{print $4}')
    job_ids+=("$job_id")

    cd ..
done

# Create the dependency string
dependency_string=$(IFS=:; echo "${job_ids[*]}")

# Create the concatenation + cleanup job script
concat_cleanup_job="concatenate_${case_name}.sh"

cat <<EOF > "$concat_cleanup_job"
#!/bin/bash
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -t 30
#SBATCH -p sapphire
#SBATCH --mem=5000
#SBATCH -o concatenate_all.out
#SBATCH -e concatenate_all.err

# Run concatenate_sorting.py
python3 $concat_script "$num_params"

# Check if the script ran successfully
if [ \$? -eq 0 ]; then
    echo "Concatenation successful, deleting num_* directories..."
    rm -rf num_*
else
    echo "Concatenation failed, directories not deleted."
    exit 1
fi
EOF

# Submit the concatenation & cleanup job only after all simulations complete successfully
sbatch --dependency=afterok:$dependency_string "$concat_cleanup_job"

echo "All jobs submitted! concatenate_all.py will run after all simulations, followed by cleanup."
