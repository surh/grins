#!/bin/bash
#SBATCH --output=Test_align.%j.out
#SBATCH --error=Test_align.%j.err
#SBATCH -p normal,hns
#SBATCH --time=500
#SBATCH --nodes=1
#SBATCH -c 6

#!/bin/bash
module load python
module load biology py-biopython
module load math py-scipystack
python ../Scripts/PKS_alignment_visualization_1x1.py -infile1 ../Test/Test.gb -infile2 ../Test/Test.gb -outfolder ../Test/ -split1 1500 -split2 1500 &
python ../Scripts/PKS_alignment_visualization_1x2.py -infile1 ../Test/Test.gb -infile2 ../Test/Test.gb -outfolder ../Test/ -split1 1500 -split2 1500 &
python ../Scripts/PKS_alignment_visualization_1x3.py -infile1 ../Test/Test.gb -infile2 ../Test/Test.gb -outfolder ../Test/ -split1 1500 -split2 1500 &
python ../Scripts/PKS_alignment_visualization_2x2.py -infile1 ../Test/Test.gb -infile2 ../Test/Test.gb -outfolder ../Test/ -split1 1500 -split2 1500 &
python ../Scripts/PKS_alignment_visualization_2x3.py -infile1 ../Test/Test.gb -infile2 ../Test/Test.gb -outfolder ../Test/ -split1 1500 -split2 1500 &
python ../Scripts/PKS_alignment_visualization_3x3.py -infile1 ../Test/Test.gb -infile2 ../Test/Test.gb -outfolder ../Test/ -split1 1500 -split2 1500 &

wait

exit 0
