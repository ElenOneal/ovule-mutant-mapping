#!/bin/bash
set -euo pipefail

# Create shell scripts to remove duplicate reads using rmdup_molbarcodes.py file


# Check if the correct number of arguments is provided
if [ "$#" -lt 3 ] || [ "$#" -gt 5 ]; then
    echo "Usage: $0 library_filenames DIR partition [time] [mem]"
    exit 1
fi

# Assign arguments to variables
file_name="$1"
DIR="$2"
partition="$3"

is_time_arg() {
  [[ "$1" =~ ^[0-9]+:[0-9]{2}:[0-9]{2}$ ]]
}

is_mem_arg() {
  [[ "$1" =~ ^[0-9]+[KMGTP]$ ]]
}

# Optional args:
# 4th arg can be either time (e.g. 24:00:00) or memory (e.g. 48G)
# 5th arg, if present, must be memory
time_limit="24:00:00"
mem_limit="4G"
if [ "$#" -ge 4 ]; then
  if is_time_arg "$4"; then
    time_limit="$4"
  elif is_mem_arg "$4"; then
    mem_limit="$4"
  fi
fi
if [ "$#" -eq 5 ]; then
  mem_limit="$5"
fi

# Check if the input file exists
if [ ! -f "$file_name" ]; then
    echo "Error: File '$file_name' not found."
    exit 2
fi

# Check if the output directory exists
if [ ! -d "$DIR" ]; then
    echo "Error: Directory '$DIR' not found."
    exit 2
fi

if ! command -v conda >/dev/null 2>&1; then
    echo "Error: conda is not available in PATH."
    exit 2
fi

script_dir="$(cd "$(dirname "$0")" && pwd)"
rmdup_script="$script_dir/rmdup_molbarcodes.py"
if [ ! -f "$rmdup_script" ]; then
  echo "Error: Helper script '$rmdup_script' not found."
  exit 2
fi

# Resolve conda base once to keep script generation fast on large sample sheets.
conda_base="$(conda info --base)"
if [ ! -f "$conda_base/etc/profile.d/conda.sh" ]; then
    echo "Error: conda initialization script not found at '$conda_base/etc/profile.d/conda.sh'."
    exit 2
fi
if ! conda env list | awk '{print $1}' | grep -Fxq 'py2'; then
    echo "Error: conda environment 'py2' not found."
    exit 2
fi

# Main loop
while IFS=$'\t' read -r a; do
  if [ -z "$a" ] || [[ "$a" =~ ^# ]]; then
    continue
  fi
  script_file="$a.rmdup.sh"
  {
    echo '#!/bin/bash'
    echo 'set -euo pipefail'
    echo '#'
    echo '#$ -S /bin/bash' 
    echo "#SBATCH --get-user-env"
    echo "#SBATCH --job-name=$a"
    echo "#SBATCH --output=${a}.rmdup.out"
    echo "#SBATCH --error=${a}.rmdup.err"
    echo '#SBATCH --cpus-per-task=1'
    echo "#SBATCH -p $partition"
    echo "#SBATCH --chdir=$DIR"
    echo "#SBATCH --mem=$mem_limit"
    echo "#SBATCH --time=$time_limit"
    echo ''
    echo "source $conda_base/etc/profile.d/conda.sh"
    echo "conda activate py2"
    echo "python \"$rmdup_script\" -p \"$a\" -s fq.gz"
  } > "$script_file"
  chmod +x "$script_file"
done < "$file_name"
