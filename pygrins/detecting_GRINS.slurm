#!/bin/bash
#SBATCH --output=Test_extract_CDS_fasta.%j.out
#SBATCH --error=Test_extract_CDS_fasta.%j.err
#SBATCH -p normal,hns
#SBATCH --time=100
#SBATCH --nodes=1
#SBATCH -c 3

#!/bin/bash

module load python
module load biology py-biopython
module load math py-scipystack
module load math py-scikit-learn

python ./detecting_GRINS_in_orphans.py -startID 0 -endID 1200 -inputf /scratch/groups/khosla/Orphan_PKS/New_approach/CDS_gb_for_GRINS/ -outputf /scratch/groups/khosla/Orphan_PKS/10AprGRINS/ &

python ./detecting_GRINS_in_orphans.py -startID 1200 -endID 2400 -inputf /scratch/groups/khosla/Orphan_PKS/New_approach/CDS_gb_for_GRINS/ -outputf /scratch/groups/khosla/Orphan_PKS/10AprGRINS/ &

python ./detecting_GRINS_in_orphans.py -startID 2400 -endID 3552 -inputf /scratch/groups/khosla/Orphan_PKS/New_approach/CDS_gb_for_GRINS/ -outputf /scratch/groups/khosla/Orphan_PKS/10AprGRINS/ &

wait

exit 0