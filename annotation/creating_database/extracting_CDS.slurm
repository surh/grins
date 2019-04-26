#!/bin/bash
#SBATCH --output=extract.%j.out
#SBATCH --error=extract.%j.err
#SBATCH -p normal,hns
#SBATCH --time=150
#SBATCH --nodes=1
#SBATCH -c 1

#!/bin/bash

module load python
module load biology py-biopython

python extracting_CDS.py -infolder /home/groups/khosla/Orphan_PKS/6dbs_prescreen/ -known ./known_clusters_MIBig_and_clusterTable.txt -outfolder ./Known_cluster_CDS/
