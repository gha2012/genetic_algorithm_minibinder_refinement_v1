#!/bin/bash
#SBATCH --gpus=a5000:1

job_id=$SLURM_JOB_ID

run_uuid=$(uuidgen)
echo "UUID for this run: $run_uuid"

file="./check.point"
if [ -f "$file" ]; then
 echo "check.point file exists. Removing it..."
 rm "$file"
else
 echo "check.point file does not exist."
fi

#This is the total number of designs that shall be created
poolFolder="pool"
iterations=50
#the lower, the more diffusion from the input structure takes place
poolSize=500
##mkdir ${outputFolder}
#change to outputFolder so that output files are not mixed up between different runs.
cd ${poolFolder}
# Set the path to the folder containing the inputfiles
input_folder_path="../input"
#put all the input .pdb files into the pool
cp $input_folder_path/*.pdb .

# Get a list of files in the folder
files=("$input_folder_path"/*.pdb)

# Check if there are any files in the folder
if [ ${#files[@]} -eq 0 ]; then
    echo "No files found in the specified folder."
    exit 1
fi

source /home/mhoe/miniconda3/etc/profile.d/conda.sh
for (( iteration=1; iteration<=$iterations; iteration++ )); do
    echo "##############inference_iteration_${iteration}#################"
    conda activate SE3nv
    # Get the number of files in the input folder
    files=(./*.pdb)
    num_files=${#files[@]}
    designsToMake=$((poolSize - num_files))
    #iterations=10
    #designsPerIteration=$(echo "scale=0; $designsToMake / $iterations" | bc)

    for (( design=1; design<=$designsToMake; design++ )); do
        echo "Iteration: $i"
        ##select a random input file for the partial diffusion
        # Generate a random index within the range of the number of files
        random_index=$((RANDOM % num_files))
        #select a random value for partial_T limited by minPartialT
        #partialT=$((RANDOM % (50 - $minPartialT + 1) + $minPartialT))
        partialT=$((RANDOM % 30 + 2))
        # Get the random file name using the random index
        random_file=${files[random_index]}
        # Print the random file name for verification (optional)
        echo "Random file selected: $random_file"
        # Store the random file name in a variable for later use
        selected_file=$random_file
        # Now you can use $selected_file in your script as needed
        #echo "Selected file stored in variable: $selected_file"
        # Count residues for chainA
        residues_chainA=$(grep '^ATOM' "$selected_file" | grep ' A ' | awk '{print $6}' | sort -u | wc -l)
        # Count residues for chainB
        residues_chainB=$(grep '^ATOM' "$selected_file" | grep ' B ' | awk '{print $6}' | sort -u | wc -l)
        first_residue_chainB=$(grep '^ATOM' "$selected_file" | grep ' B ' | awk 'NR==1{print $6}')
        last_residue_chainB=$(grep '^ATOM' "$selected_file" | grep ' B ' | awk 'END{print $6}')
        contigString="[${residues_chainA}-${residues_chainA}/0 B${first_residue_chainB}-${last_residue_chainB}]"
        /home/${USER}/RFdiffusion/scripts/run_inference.py inference.output_prefix="./iteration_${iteration}_design_${design}" inference.input_pdb="./${selected_file}" contigmap.contigs="${contigString}" inference.num_designs=1 diffuser.partial_T="${partialT}"
    done

    conda deactivate
    conda activate proteinmpnn_binder_design
    echo "##############convert_iteration_${iteration}###################"
    silentfrompdbs iteration_${iteration}_design_*.pdb > "inference_output_iteration_${iteration}.silent"
    #now that all has been stored in a silent file, delete the new pdb files, which only contain poly Gly models.
    find . -maxdepth 1 -type f -name "iteration_${iteration}_design_*.pdb" ! -name "$inputPdb" -exec rm {} \;
    rm *.trb

    echo "############interface design_iteration_${iteration}############"
    /home/${USER}/dl_binder_design/mpnn_fr/dl_interface_design.py -silent "inference_output_iteration_${iteration}.silent"
    mv out.silent "interfaceDesignOutput_iteration_${iteration}.silent"

    echo "############alphaFold check_iteration_${iteration}############"
    conda deactivate
    conda activate af2_binder_design
    /home/${USER}/dl_binder_design/af2_initial_guess/predict.py -silent "interfaceDesignOutput_iteration_${iteration}.silent"
    mv out.silent "afCheck_iteration_${iteration}.silent"
    mv out.sc "afCheck_iteration_${iteration}.sc"

    echo "###########remove bad designs_iteration_${iteration}##########"
    #extract the scored files
    silentextract "./afCheck_iteration_${iteration}.silent"
    #keep only top 50% of designs
    declare -a files_array
    # Iterate through each PDB file in the pool and extract the pae_interaction values
    for file in *.pdb; do
        pae_value=$(grep 'pae_interaction' "$file" | awk '{print $NF}')
        files_array+=("$file $pae_value")
    done

    # Sort the array based on pae_interaction values (ascending order)
    IFS=$'\n' sorted_array=($(sort -k2 -n <<<"${files_array[*]}"))
    unset IFS

    # Calculate the index for the 50% worst files
    index=$(((${#sorted_array[@]} * 50) / 100))
    lengthOfArray=${#sorted_array[@]}
    echo $index

    # Iterate through the 50% worst files and delete them
    for ((i = index; i < $lengthOfArray; i++)); do
        file_to_delete=$(echo "${sorted_array[i]}" | awk '{print $1}')
        rm "$file_to_delete"
        echo "Deleted: $file_to_delete"
    done
    #now only the "good" files remain. go for the next iteration.
    conda deactivate
done

#echo "############cleanup##################"
conda deactivate
mv "../slurm-${job_id}.out" .
