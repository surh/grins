import os
import argparse

def process_arguments():
    # Read arguments
    parser_format = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(formatter_class=parser_format)
    required = parser.add_argument_group("Required arguments")


    # Define required arguments
    required.add_argument("--input_file", help=("File with the results of genome clustering."),
                          required=True, type=str)

    required.add_argument("--output_folder", help=("Folder to store files for further GRINS detection"),
                          required=True, type=str)


    args = parser.parse_args()

    return args


args = process_arguments()

with open(args.input_file,'r') as genome_clustering_results_file:
	genome_clustering_results=genome_clustering_results_file.readlines()

genome_clusters={}
for i in range(1,len(genome_clustering_results)):
	line_data=genome_clustering_results[i].strip("\n").split("\t")
	genome=line_data[0]
	cluster=line_data[1]
	if cluster in genome_clusters:
		genome_clusters[cluster]=genome_clusters[cluster].append(genome)
	else:
		genome_clusters[cluster]=[genome]

chosen_genomes=[]
for cluster_key in genome_clusters:
	curr_clusters=genome_clusters[cluster_key]
	curr_genome=np.random.choice(curr_clusters)
	# curr_genome_renamed=curr_genome[:curr_genome.find(".f")]+".fasta"
	chosen_genomes.append(curr_genome_renamed)

os.system("cp "+" ".join(chosen_genomes)+" "+args.output_folder)