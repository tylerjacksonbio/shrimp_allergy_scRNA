#!/bin/sh
#SBATCH --mail-type=end
#SBATCH --mail-type=fail
#SBATCH --mail-user=tyler.jackson@bcm.edu
#SBATCH --job-name=cellranger_batch
#SBATCH --nodes=1
#SBATCH --cpus-per-task=25
#SBATCH --time=70:00:00
#SBATCH --mem=96GB
#SBATCH -e /mount/hli/Tyler/scRNA_Scripts/scRNA_error_files/%j.e
#SBATCH -o /mount/hli/Tyler/scRNA_Scripts/scRNA_output_files/%j.o

### Commonly used parameters

### Settings of commonly used folders
mountFolder=/mount/hli
projectFolder=/project/hli
codeFolder=/home/u239500/codeFolder/HuiLab

############### Main Code ###############

# Load cellranger
module load cellranger/6.1.2

# Setting folder paths
analysisFolder=/mount/hli/Tyler/Shrimp_Allergy_Collaboration/Cellranger_Outputs
# Folder of the indexed genome
index10xFolder=/mount/hli/Tyler/Homo_sapiens_reference/refdata_gex_GRCh38_2020_A

# Define the lists
sampleNames=("CombinedESC24" "CombinedESC25" "CombinedESC26" "CombinedESC27")
folderNames=("SD1" "SD2" "SD3" "SD4")

### Loop through the lists
for i in ${!sampleNames[@]}; do
    sampleName=${sampleNames[$i]}
    folderName=${folderNames[$i]}
    cellrangerFolder=$analysisFolder/$folderName

    # Create folder if it doesn't exist
    mkdir -p $cellrangerFolder

    # Set fastq folder
    fastqFolder=$mountFolder/Tyler/Shrimp_Allergy_Collaboration/fastq_files

    # Change directory to the folder
    cd $cellrangerFolder

    # Run cellranger count
    cellranger count --id=$folderName \
        --transcriptome=$index10xFolder \
        --chemistry=ARC-v1 \
        --fastqs=$fastqFolder \
        --sample=$sampleName \
        --localmem=96 \
        --include-introns=true \
        --create-bam=true
done

# Unload cellranger module
module unload cellranger
