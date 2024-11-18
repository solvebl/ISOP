# READS TO TREES WORKFLOW v. 160223

This workflow is meant for handling target capture data using the angiosperms353 bait kit. The input data is raw reads (.fastq) and final output is a species tree from ASTRAL (.tree). The workflow should be run on the Saga cluster (https://documentation.sigma2.no/hpc_machines/saga.html). All code starting with #!/bin/bash should be saved as separate, plain text files, and run as a jobscripts (https://documentation.sigma2.no/jobs/job_scripts.html).

# Before you start

## Log in to Saga 
```
ssh [your_username]@saga.sigma2.no
```
## Copy template directory
```
cd /cluster/projects/nn9835k
cp -r msc_template [your_name]
```

## Move raw files to your user directory 
```
mv [full_path_to_raw_files_on_saga]/*.fastq) [your_name]/raw_reads
cd [your_name] 
```

## Edit your namelist
The namelist should contain a list of the unique file/sample prefixes of your raw files, one per line. Also, edit raw file names if necessary.

```
vi jobscripts/namelist 
```

Now, the following steps can be run from your user directory (cluster/projects/nn9835k/[your_name]. N is always the number of samples you have, and the number of rows in namelist.


# 1 - Unzip fastq files

```
sbatch --array=1-N jobscripts/unzip.sh
```
```
#!/bin/bash
#SBATCH --account=nn9835k
#SBATCH --job-name=unzip
#SBATCH --time=06:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --ntasks=1 --cpus-per-task=1 --ntasks-per-node=1

#set name variable
CURRENT_NAME=$(awk NR==$SLURM_ARRAY_TASK_ID jobscripts/namelist)
R1=raw_reads/$CURRENT_NAME"_R1.fastq.gz"
R2=raw_reads/$CURRENT_NAME"_R2.fastq.gz"

#execute
gunzip -r $R1
gunzip -r $R2
```

# 2 - FastQC on raw reads 
Check quality of raw reads using the program FastQC (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/). 

```
sbatch --array=1-N jobscripts/fastqc_raw.sh
```
```
#!/bin/bash
#SBATCH --account=nn9835k
#SBATCH --job-name=fastqc_raw
#SBATCH --time=06:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --ntasks=1 --cpus-per-task=1 --ntasks-per-node=1

#clear any inherited modules
module --quiet purge 
#exit on errors
set -o errexit 
set -o nounset

#load modules
module load Python/3.8.2-GCCcore-9.3.0
module load FastQC/0.11.9-Java-11
module list
	
#set name variable
CURRENT_NAME=$(awk NR==$SLURM_ARRAY_TASK_ID jobscripts/namelist)
R1=raw_reads/$CURRENT_NAME"_R1.fastq"
R2=raw_reads/$CURRENT_NAME"_R2.fastq"

#execute
fastqc -o results_quality_raw  $R1
fastqc -o results_quality_raw  $R2
```
Transfer files to computer for visual check. This command should not be run from Saga, but from a new Ubuntu window (or terminal window on Mac).
```
#ubuntu
scp -r [your_username]@saga.sigma2.no:[full_path_on_saga]/results_quality_raw/*.html /mnt/c/ubuntu/[your_directory]
#mac
scp -r [your_username]@saga.sigma2.no:[full_path_on_saga]/results_quality_raw/*.html [path_on_mac]
```

# 3 - Clean reads
Use the program Trimmomatic (http://www.usadellab.org/cms/?page=trimmomatic) to remove adapters, low quality bases etc. Read manual (http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf) for detailed descriptions of the different settings.

```
sbatch --array=1-N jobscripts/trimmomatic.sh 
```
```
#!/bin/bash
#SBATCH --account=nn9835k
#SBATCH --job-name=trimmomatic
#SBATCH --time=16:00:00
#SBATCH --mem-per-cpu=12G
#SBATCH --ntasks=1 --cpus-per-task=1 --ntasks-per-node=1

#clear any inherited modules
module --quiet purge 
module load Java/17.0.2
# exit on errors
set -o errexit
set -o nounset

#set name variable
CURRENT_NAME=$(awk NR==$SLURM_ARRAY_TASK_ID jobscripts/namelist)
R1=raw_reads/$CURRENT_NAME"_R1.fastq"
R2=raw_reads/$CURRENT_NAME"_R2.fastq"

#execute
java -jar /cluster/software/Trimmomatic/0.39-Java-11/trimmomatic.jar PE "$R1" "$R2" cleaned_reads/"$CURRENT_NAME"_R1_paired.fastq  cleaned_reads/"$CURRENT_NAME"_R1_unpaired.fastq  cleaned_reads/"$CURRENT_NAME"_R2_paired.fastq cleaned_reads/"$CURRENT_NAME"_R2_unpaired.fastq  ILLUMINACLIP:/cluster/software/Trimmomatic/0.39-Java-11/adapters/TruSeq3-PE-2.fa:2:30:10  LEADING:30 TRAILING:30 SLIDINGWINDOW:4:30 MINLEN:36

#zip raw reads after cleaning
gzip $R1
gzip $R2
```

# 4 - FastQC on cleaned reads	
Check quality of cleaned reads using the program FastQC (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).

```
sbatch --array=1-N jobscripts/fastqc_cleaned.sh
```
```
#!/bin/bash
#SBATCH --account=nn9835k
#SBATCH --job-name=fastqc_clean
#SBATCH --time=06:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --ntasks=1 --cpus-per-task=1 --ntasks-per-node=1

#clear any inherited modules
module --quiet purge 
# exit on errors
set -o errexit 
set -o nounset

#load modules
module load Python/3.8.2-GCCcore-9.3.0
module load FastQC/0.11.9-Java-11
module list
	
#set name variable
CURRENT_NAME=$(awk NR==$SLURM_ARRAY_TASK_ID jobscripts/namelist)
R1=cleaned_reads/$CURRENT_NAME"_R1_paired.fastq"
R2=cleaned_reads/$CURRENT_NAME"_R2_paired.fastq"

#execute
fastqc -o results_quality_cleaned  $R1
fastqc -o results_quality_cleaned  $R2
```
Transfer files to computer for visual check. This command should not be run from Saga, but from a new Ubuntu window (or terminal window on Mac).
```
#ubuntu
scp -r [your_username]@saga.sigma2.no:[full_path_on_saga]/results_quality_cleaned/*.html /mnt/c/ubuntu/[your_directory]
#mac
scp -r [your_username]@saga.sigma2.no:[full_path_on_saga]/results_quality_cleaned/*.html [path_on_mac]
```
If quality is good, move on to next step. If not, re-run step 3 with different settings (see Trimmomatic manual http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf).

# 5 - Zip cleaned reads

```
sbatch --array=1-N jobscripts/zip_cleaned.sh
```
```
#!/bin/bash
#SBATCH --account=nn9835k
#SBATCH --job-name=zip
#SBATCH --time=06:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --ntasks=1 --cpus-per-task=1 --ntasks-per-node=1

#set name variable
CURRENT_NAME=$(awk NR==$SLURM_ARRAY_TASK_ID jobscripts/namelist)
R1=cleaned_reads/$CURRENT_NAME"_R1_paired.fastq"
R2=cleaned_reads/$CURRENT_NAME"_R2_paired.fastq"
R3=cleaned_reads/$CURRENT_NAME"_R1_unpaired.fastq"
R4=cleaned_reads/$CURRENT_NAME"_R2_unpaired.fastq"

#execute
gzip $R1
gzip $R2
gzip $R3
gzip $R4
```

# 6 - HybPiper assemble
Use the pipeline Hybpiper to assemble and map your clean reads to a reference (https://github.com/mossmatters/HybPiper/wiki/Tutorial). This part of the pipeline generates very many temporary and output files. If you have many samples, you may have to submit your samples in bulk (e.g. --array=1-35, --array=36-70 etc.), to avoid exceeding disk quota. Wait for one bulk to finish before you submit a new one. Extraction of introns is optional (https://github.com/mossmatters/HybPiper/wiki/Introns).

```
sbatch --array=1-N jobscripts/hybpiper_assemble.sh
```
```
#!/bin/bash
#SBATCH --account=nn9835k 
#SBATCH --time=12:00:00
#SBATCH --job-name=hybpiper_assemble
#SBATCH --nodes=1 --ntasks=8
#SBATCH --mem-per-cpu=4G

#load Anaconda3
module load Anaconda3/2019.03

#set the ${PS1} (needed in the source of the Anaconda environment)
export PS1=\$

#source the conda environment setup
source ${EBROOTANACONDA3}/etc/profile.d/conda.sh
#source ${EBROOTMINICONDA3}/etc/profile.d/conda.sh

#deactivate any spill-over environment from the login node
conda deactivate &>/dev/null

#activate the environment by using the full path (not name)
conda activate /cluster/projects/nn9835k/conda_envs/hybpiper
cd hybpiper

#set name variable
CURRENT_NAME=$(awk NR==$SLURM_ARRAY_TASK_ID ../jobscripts/namelist)
R1=../cleaned_reads/$CURRENT_NAME"_R1_paired.fastq.gz"
R2=../cleaned_reads/$CURRENT_NAME"_R2_paired.fastq.gz"

#execute
hybpiper assemble -r $R1 $R2 -t_dna /cluster/projects/nn9835k/Solveig/NewTargets/mega353.fasta --bwa --prefix $CURRENT_NAME --cpu=${SLURM_NTASKS}  
#--run_intronerate
#--run_intronerate --start_from exonerate_contigs
```

# 7 - HybPiper stats
Get statistics for the HybPiper run. 

```
sbatch jobscripts/hybpiper_stats.sh
```
```
#!/bin/bash
#SBATCH --account=nn9835k 
#SBATCH --time=6:00:00
#SBATCH --job-name=hybpiper_stats
#SBATCH --nodes=1 --ntasks=8
#SBATCH --mem-per-cpu=4G

#load Anaconda3
module load Anaconda3/2019.03

#set the ${PS1} (needed in the source of the Anaconda environment)
export PS1=\$

#source the conda environment setup
source ${EBROOTANACONDA3}/etc/profile.d/conda.sh

#deactivate any spill-over environment from the login node
conda deactivate &>/dev/null

#activate the environment by using the full path (not name)
conda activate /cluster/projects/nn9835k/conda_envs/hybpiper
cd hybpiper 

#execute
hybpiper stats -t_dna /cluster/projects/nn9835k/Solveig/NewTargets/mega353.fasta gene ../jobscripts/namelist
```

# 8 - HybPiper recovery heatmap 
Make heatmap to summarize HybPiper run.

```
sbatch jobscripts/hybpiper_recovery_heatmap.sh
```
```
#!/bin/bash
#SBATCH --account=nn9835k 
#SBATCH --time=6:00:00
#SBATCH --job-name=hybpiper_stats
#SBATCH --nodes=1 --ntasks=8
#SBATCH --mem-per-cpu=4G

#load Anaconda3
module load Anaconda3/2019.03

#set the ${PS1} (needed in the source of the Anaconda environment)
export PS1=\$

#source the conda environment setup
source ${EBROOTANACONDA3}/etc/profile.d/conda.sh

#deactivate any spill-over environment from the login node
conda deactivate &>/dev/null

#activate the environment by using the full path (not name)
conda activate /cluster/projects/nn9835k/conda_envs/hybpiper
cd hybpiper

#execute
hybpiper recovery_heatmap seq_lengths.tsv
```

# 9 - HybPiper retrieve sequences
Retrieve the unaligned sequences as fasta files, one file per gene.

```
sbatch jobscripts/hybpiper_retrieve_sequences.sh
```
```
#!/bin/bash
#SBATCH --account=nn9835k 
#SBATCH --time=6:00:00
#SBATCH --job-name=hybpiper_retrieve_seqs
#SBATCH --nodes=1 --ntasks=8
#SBATCH --mem-per-cpu=4G

#load Anaconda3
module load Anaconda3/2019.03

#set the ${PS1} (needed in the source of the Anaconda environment)
export PS1=\$

#source the conda environment setup
source ${EBROOTANACONDA3}/etc/profile.d/conda.sh

#deactivate any spill-over environment from the login node
conda deactivate &>/dev/null

#activate the environment by using the full path (not name)
conda activate /cluster/projects/nn9835k/conda_envs/hybpiper
cd hybpiper

#execute
hybpiper retrieve_sequences dna -t_dna /cluster/projects/nn9835k/Solveig/NewTargets/mega353.fasta --sample_names ../jobscripts/namelist --fasta_dir sequences
```

The next two steps are optional.

# 10 - Retrieve paralogs
Retrieve paralogous sequences (https://github.com/mossmatters/HybPiper/wiki/Paralogs)

```
sbatch jobscripts/hybpiper_paralogs.sh
```
```
#!/bin/bash
#SBATCH --account=nn9835k 
#SBATCH --time=6:00:00
#SBATCH --job-name=hybpiper_paralogs
#SBATCH --nodes=1 --ntasks=8
#SBATCH --mem-per-cpu=4G

#load Anaconda3
module load Anaconda3/2019.03

#set the ${PS1} (needed in the source of the Anaconda environment)
export PS1=\$

#source the conda environment setup
source ${EBROOTANACONDA3}/etc/profile.d/conda.sh

#deactivate any spill-over environment from the login node
conda deactivate &>/dev/null

#activate the environment by using the full path (not name)
conda activate /cluster/projects/nn9835k/conda_envs/hybpiper
cd hybpiper

#execute 
hybpiper paralog_retriever ../jobscripts/namelist -t_dna /cluster/projects/nn9835k/Solveig/NewTargets/mega353.fasta
```
View the output file 'paralogs_above_threshold_report.txt', select and copy the list of genes. Open a plain text file called 'paraloglist'. Paste the list of genes and save the file. Move the new file to your HybPiper directory.
```
cat paralogs_above_threshold_report.txt
vi paraloglist
mv paraloglist hybpiper
```

# 11 - Make paralog trees
Make gene trees of all paralogs to view the relatedness of copies.
Example 1: two copies from the same species form a clade --> no sign of ancient duplication event.
Example 2: a sequence is more related to a copy in a different species than it is to a copy in the same species --> sign of ancient duplication event.
Here, N is the number of genes with paralog warnings, and the number of rows in paraloglist.

```
sbatch --array=1-N jobscripts/paralog_tree.sh
```
```
#!/bin/bash
#SBATCH --account=nn9835k
#SBATCH --job-name=paralog_gene_tree
#SBATCH --time=06:00:00
#SBATCH --mem-per-cpu=12G
#SBATCH --ntasks=1 --cpus-per-task=1 --ntasks-per-node=1

#clear any inherited modules
module --quiet purge  
#load modules
module load MAFFT/7.490-GCC-10.3.0-with-extensions
module load FastTree/2.1.11-GCCcore-10.3.0 
module list

cd hybpiper

#set the input file to process
name=$(sed -n ${SLURM_ARRAY_TASK_ID}p paraloglist)

#execute 
cat paralogs_all/"$name"_paralogs_all.fasta | mafft --auto - | FastTree -nt -gtr > "$name".paralogs.tre
```
To clean up, create a paralog directory, and move all paralog-files there.
```
mkdir paralogfiles
mv *.tre paralogfiles
mv paralogs* paralogfiles
mv paralog_* paralogfiles
mv paraloglist paralogfiles
```
Move whatever you need to computer. Run this command from a new Ubuntu window (or terminal window on Mac).
```
#ubuntu
scp -r [your_username]@saga.sigma2.no:/cluster/projects/nn9835k/[your_name]/hybpiper /mnt/c/ubuntu/[your_directory]
#mac
scp -r [your_username]@saga.sigma2.no:/cluster/projects/nn9835k/[your_name]/hybpiper [path_on_mac]
```
Inspect your parlog trees manually. For all genes that show signs of an ancient duplication event, remove them form hybpiper/sequences/ before continuing. 

# 12 - Download outgroups
Download outgroups from kew tree of life (https://treeoflife.kew.org/tree-of-life).
Edit the jobscript according to which genera you want to include. 
Prepare a genelist the same way as you prepared your namelist: a list of genes names, one per line. Make sure to edit your genelist if you removed some genes from the analysis following paralog investigation in step 11.

```
sbatch --array=1-353 jobscripts/get_outgroups.sh 
```
```
#!/bin/bash
#SBATCH --account=nn9835k
#SBATCH --job-name=get_outgroups
#SBATCH --time=01:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --ntasks=1 --cpus-per-task=1 --ntasks-per-node=1

#load Anaconda3
module load Anaconda3/2019.03
module load Biopython/1.79-foss-2022a

# set the input file to process
name=$(sed -n ${SLURM_ARRAY_TASK_ID}p jobscripts/genelist)

#execute
python /cluster/projects/nn9835k/kew_353_outgroups/breakfasta.py -f /cluster/projects/nn9835k/kew_353_outgroups/"$name".dna.fasta -d "$name"_outgroups.fasta -i Brodiaea Ornithogalum -s

mv "$name"_outgroups.fasta outgroups
```

# 13 - Concatenate in- and outgroup
Concatenate outgroup-fasta-files with your ingroup-fasta-files. 

```
sbatch --array=1-353 jobscripts/concat_fastas.sh 
```
```
#!/bin/bash
#SBATCH --account=nn9835k
#SBATCH --job-name=concat_fastas
#SBATCH --time=01:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --ntasks=1 --cpus-per-task=1 --ntasks-per-node=1

#load Anaconda3
module load Anaconda3/2019.03
module load Biopython/1.79-foss-2022a

#set the input file to process
name=$(sed -n ${SLURM_ARRAY_TASK_ID}p jobscripts/genelist)

#execute
python3 /cluster/projects/nn9835k/conda_envs/AMAS/lib/python3.10/site-packages/amas/AMAS.py concat -f fasta -d dna -i hybpiper/sequences/"$name".FNA outgroups/"$name"_outgroups.fasta -u fasta -t "$name"_w_outgroups.fasta

mv "$name"_w_outgroups.fasta unaligned
```

# 14 - Change header names inside fasta-file
We need to change the naming of the headers inside the fasta files for the rest of the steps to work properly.
Change input and output folders in python script before running.

```
python3 change_fastas.py
```
```
import os

#define the directory containing the input files
input_dir = '/cluster/projects/nn9835k/Solveig/unaligned_fastas/kew_outgroups'

#define the directory to write the updated files to
output_dir = '/cluster/projects/nn9835k/Solveig/unaligned_fastas/kew_outgroups/edited_fastas'

#loop over each file in the input directory
for filename in os.listdir(input_dir):
    #only consider files with the .fasta extension
    if filename.endswith(".fasta"):
        #open the input file for reading and the output file for writing
        with open(os.path.join(input_dir, filename), "r") as input_file, \
             open(os.path.join(output_dir, filename), "w") as output_file:
            #loop over each line in the input file
            for line in input_file:
                #check if the line is a header line
                if line.startswith(">"):
                    #split the header line by whitespace
                    header_parts = line.strip().split()
                    #only keep the first part of the header
                    header = header_parts[0]
                    #write out the updated header to the output file
                    output_file.write(header + "\n")
                else:
                    #write out the sequence line to the output file
                    output_file.write(line)
```

# 15 - Multiple alignments and alignment trimming 
Alignments are made with mafft (https://mafft.cbrc.jp/alignment/software/linux.html), trimming with clipkit (https://github.com/JLSteenwyk/ClipKIT).
Settings might need adjustment, check software documentation for details.

```
sbatch --array=1-353 jobscripts/align_trim.sh
```
```
#!/bin/bash
#SBATCH --account=nn9835k 
#SBATCH --time=6:00:00
#SBATCH --job-name=align_trim
#SBATCH --nodes=1 --ntasks=8
#SBATCH --mem-per-cpu=4G

#load Anaconda3
module load Anaconda3/2019.03

#set the ${PS1} (needed in the source of the Anaconda environment)
export PS1=\$

#source the conda environment setup
source ${EBROOTANACONDA3}/etc/profile.d/conda.sh
#source ${EBROOTMINICONDA3}/etc/profile.d/conda.sh

#deactivate any spill-over environment from the login node
conda deactivate &>/dev/null

conda activate /cluster/projects/nn9835k/conda_envs/clipkit
module --quiet purge
module load MAFFT/7.470-gompi-2020a-with-extensions

#set the input file to process
name=$(sed -n ${SLURM_ARRAY_TASK_ID}p jobscripts/genelist)

#run mafft 
mafft --auto unaligned/"$name"_w_outgroups_edited.fasta > "$name"_aligned.fasta
mv "$name"_aligned.fasta aligned
#run clipcit
clipkit 
clipkit aligned/"$name"_aligned.fasta -o "$name"_aligned_trimmed.fasta -m gappy
mv "$name"_aligned_trimmed.fasta aligned_trimmed
```

# 16 - Alignment summaries
Summarizes alignment statistics using AMAS (https://github.com/marekborowiec/AMAS).

```
sbatch jobscripts/alignment_summary.sh
```
```
#!/bin/bash
#SBATCH --account=nn9835k
#SBATCH --job-name=alignment_summary
#SBATCH --time=01:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --ntasks=1 --cpus-per-task=1 --ntasks-per-node=1

#load Anaconda3
module load Anaconda3/2019.03
module load Biopython/1.79-foss-2022a

#execute
python3 /cluster/projects/nn9835k/conda_envs/AMAS/lib/python3.10/site-packages/amas/AMAS.py summary -f fasta -d dna -i aligned/*.fasta -o alignment_summaries.txt
```

# 17 - Make gene trees
Make individual gene trees for all your trimmed alignments using iq-tree (http://www.iqtree.org/).
You might want to adjust settings, see software documentation for details.

```
sbatch --array=1-353 jobscripts/iqtree.sh
```
```
#!/bin/bash
#SBATCH --account=nn9835k 
#SBATCH --time=6:00:00
#SBATCH --job-name=iqtree_genetrees
#SBATCH --nodes=1 --ntasks=8
#SBATCH --mem-per-cpu=4G

module --quiet purge
module load Anaconda3/2022.05
module load IQ-TREE/2.2.1-gompi-2021b 
module list

#set the input file to process
name=$(sed -n ${SLURM_ARRAY_TASK_ID}p jobscripts/genelist)

#execute
iqtree2 -s aligned_trimmed/"$name"_aligned_trimmed.fasta -B 1000
mv aligned_trimmed/"$name"_aligned_trimmed.fasta.* iqtree
```

# 18 - Concatenate gene trees
Concatenate all gene tree files as required for the neext step.

```
sbatch jobscripts/concat_genetrees.sh
```
```
#!/bin/bash
#SBATCH --account=nn9835k
#SBATCH --job-name=concat_gene_trees
#SBATCH --time=01:00:00
#SBATCH --mem-per-cpu=3G
#SBATCH --ntasks=1 --cpus-per-task=1 --ntasks-per-node=1

cat iqtree/*.contree > astral/genetrees_concat.trees
```

# 19 - Make species tree 
Make species tree using astral (https://github.com/smirarab/ASTRAL/blob/master/astral-tutorial.md).
You might want to adjust settings, see software documentation for details.               

```
sbatch jobscripts/astral.sh
```
```
#!/bin/bash
#SBATCH --account=nn9835k 
#SBATCH --time=8:00:00
#SBATCH --job-name=astral
#SBATCH --nodes=1 --ntasks=8
#SBATCH --mem-per-cpu=4G

#clear any inherited modules
module --quiet purge 
#load java
module load Java/17.0.2
#exit on errors
set -o errexit
set -o nounset

java -jar /cluster/projects/nn9835k/Astral/astral.5.7.8.jar -i genetrees_concat.trees -o speciestree.tree 2> astral.log
```

# 20 - Change tip labels 
Changes names on tree tips, according to names of your accessions.
You need to provide a key file (key_file.txt) with current and new names of samples/tips separated by tab, each sample on a new row. 

```
sbatch jobscripts/change_tip_labels.sh
```
```
#!/bin/bash
#SBATCH --account=nn9835k
#SBATCH --job-name=change_fastaheaders
#SBATCH --time=01:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --ntasks=1 --cpus-per-task=1 --ntasks-per-node=1

#load Anaconda3
module load Anaconda3/2019.03

#set the ${PS1} (needed in the source of the Anaconda environment)
export PS1=\$

#source the conda environment setup
source ${EBROOTANACONDA3}/etc/profile.d/conda.sh
#source ${EBROOTMINICONDA3}/etc/profile.d/conda.sh

#deactivate any spill-over environment from the login node
conda deactivate &>/dev/null
#activate ete3
conda activate /cluster/projects/nn9835k/conda_envs/ete3

#execute
python3 /cluster/projects/nn9835k/kew_353_outgroups/change_tip_labels.py astral/species.tree jobscripts/key_file.txt > astral_tree_edited.tree
mv astral_tree_edited.tree astral
```

Now, move your edited tree-file to your computer, and look at the beautiful tree you made!!
