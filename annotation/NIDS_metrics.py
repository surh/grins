import os
import numpy as np
import scipy.stats
from scipy.stats.stats import pearsonr
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.colors as clr
import seaborn as sns
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.SeqFeature import SeqFeature,FeatureLocation

sns.set(style="whitegrid",font_scale=1.5)

font = {'family': 'sans-serif',
        'color':  'black',
        'weight': 'normal',
        'size': 20, 'rotation':0,
        'verticalalignment': 'bottom',
        'horizontalalignment': 'left'
        }

GRINS_file=open("NIDS_GRINS_list.txt",'w')

# Calculating skews
def GC_skew(string):
    G_number=string.count("G")
    C_number=string.count("C")
    skew=float(G_number-C_number)/float(G_number+C_number)
    GCcontent=float(G_number+C_number)/float(len(string))
    return G_number,C_number,skew,GCcontent
def TA_skew(string):
    A_number=string.count("A")
    T_number=string.count("T")
    skew=float(T_number-A_number)/float(A_number+T_number)
    return A_number,T_number,skew
def AT_skew(string):
    A_number=string.count("A")
    T_number=string.count("T")
    skew=float(A_number-T_number)/float(A_number+T_number)
    return A_number,T_number,skew


# Calculating p-value to display on graphs
def p_value_def(p_value):
    if p_value<0.001:
        return_p_val="\np<0.001"
    elif p_value<0.01:
        return_p_val="\np<0.01"
    elif p_value<0.05:
        return_p_val="\np<0.05"
    else:
        return_p_val="\nns (p>0.05)"
    return return_p_val

# Transposing a matrix

def transpose_matrix(data):
	data_transposed=[]
	for j in range(0,len(data[0])):
		data_transposed.append([])
	for i in range(0,len(data)):
		for j in range(0,len(data[i])):
			data_transposed[j].append(data[i][j])
	return data_transposed


### NIDS alignmemt scores ###

# Getting  scores for NIDS alignment to itself
with open("./NIDS_vs_NIDS_window150_step30_1x1.txt", 'r') as input_file1:
	input_data1=input_file1.readlines()
for i in range(0,len(input_data1)):
	line=input_data1[i].split(",")
	input_data1[i]=line[:-1]      
with open("./NIDS_vs_NIDS_window150_step30_1x2.txt", 'r') as input_file2:
	input_data2=input_file2.readlines()
for i in range(0,len(input_data2)):
	line=input_data2[i].split(",")
	input_data2[i]=line[:-1]
with open("./NIDS_vs_NIDS_window150_step30_1x3.txt", 'r') as input_file3:
	input_data3=input_file3.readlines()
for i in range(0,len(input_data3)):
	line=input_data3[i].split(",")
	input_data3[i]=line[:-1]
with open("./NIDS_vs_NIDS_window150_step30_2x2.txt", 'r') as input_file4:
	input_data4=input_file4.readlines()
for i in range(0,len(input_data4)):
	line=input_data4[i].split(",")
	input_data4[i]=line[:-1]
with open("./NIDS_vs_NIDS_window150_step30_2x3.txt", 'r') as input_file5:
	input_data5=input_file5.readlines()
for i in range(0,len(input_data5)):
	line=input_data5[i].split(",")
	input_data5[i]=line[:-1]
with open("./NIDS_vs_NIDS_window150_step30_3x3.txt", 'r') as input_file6:
	input_data6=input_file6.readlines()
for i in range(0,len(input_data6)):
	line=input_data6[i].split(",")
	input_data6[i]=line[:-1]
input_data2_transposed=transpose_matrix(input_data2)
input_data3_transposed=transpose_matrix(input_data3)
input_data5_transposed=transpose_matrix(input_data5)

# Constructing the matrix of alignment scores
homology_matrix=[]
for i in range(0, len(input_data1)):
	homology_matrix.append([])
	line_data1=input_data1[i]
	line_data2=input_data2[i]
	line_data3=input_data3[i]
	for j in range(0,len(line_data1)):
		homology_matrix[len(homology_matrix)-1].append(float(line_data1[j])/1.5)
	for j in range(0,len(line_data2)):
		homology_matrix[len(homology_matrix)-1].append(float(line_data2[j])/1.5)
	for j in range(0,len(line_data3)):
		homology_matrix[len(homology_matrix)-1].append(float(line_data3[j])/1.5)
for i in range(0,len(input_data2_transposed)):
	homology_matrix.append([])
	line_data2=input_data2_transposed[i]
	line_data4=input_data4[i]
	line_data5=input_data5[i]
	for j in range(0,len(line_data2)):
		homology_matrix[len(homology_matrix)-1].append(float(line_data2[j])/1.5)
	for j in range(0,len(line_data4)):
		homology_matrix[len(homology_matrix)-1].append(float(line_data4[j])/1.5)
	for j in range(0,len(line_data5)):
		homology_matrix[len(homology_matrix)-1].append(float(line_data5[j])/1.5)
for i in range(0,len(input_data3_transposed)):
 	homology_matrix.append([])
 	line_data3=input_data3_transposed[i]
 	line_data5=input_data5_transposed[i]
 	line_data6=input_data6[i]
 	for j in range(0,len(line_data3)):
 		homology_matrix[len(homology_matrix)-1].append(float(line_data3[j])/1.5)
 	for j in range(0,len(line_data5)):
 		homology_matrix[len(homology_matrix)-1].append(float(line_data5[j])/1.5)
 	for j in range(0,len(line_data6)):
 		homology_matrix[len(homology_matrix)-1].append(float(line_data6[j])/1.5)   

# Constructing a list of maximum identity scores
similarity_max_array=[]
for j in range(0,len(homology_matrix)):
    if j>5 and j<len(homology_matrix[j])-6:
        max1=np.max(homology_matrix[j][:j-5])
        max2=np.max(homology_matrix[j][j+6:])
        similarity_max_array.append(np.max([max1,max2]))
    elif j<=5:
        similarity_max_array.append(np.max(homology_matrix[j][j+6:]))
    else:
        similarity_max_array.append(np.max(homology_matrix[j][:j-5]))

# Determining GRINS regions, based on max identity scores, and parameters
region_starts=[]
region_ends=[]
possible_region_starts=[]
possible_region_ends=[]

local_start_cutoff=75
next10_start_cutoff=85
local_end_cutoff=80
previous10_end_cutoff=80

for i in range(5,len(similarity_max_array)):
    if len(region_starts)==len(region_ends):
        if similarity_max_array[i]>=local_start_cutoff and np.mean(similarity_max_array[i:i+10])>=next10_start_cutoff:
            region_starts.append(i)
    else:
        if np.mean(similarity_max_array[i-11:i-1])>=previous10_end_cutoff and similarity_max_array[i]<local_end_cutoff:
            region_ends.append(i)

print region_starts
print region_ends


# # OPTIONAL: Manually adding GRINS regions, for example if an apparent GRINS region is present in only 1 copy
# GRINS_file.write("Manually added GRINSs: \n")
# region_starts.append(358)
# region_starts.append(683)
# region_ends.append(400)
# region_ends.append(720)
# GRINS_file.write(str(358*30)+"-"+str(400*30)+"\n")
# GRINS_file.write(str(683*30)+"-"+str(720*30)+"\n")
# GRINS_file.write("Manually added erased/possible GRINSs (GRINS?)\n")
# possible_region_starts=[400]
# possible_region_ends=[434]
# region_starts=[]
# region_ends=[]

### NIDS skew values ###

# Reading the GenBank file
record = SeqIO.read("NIDS.gb", "gb")
input_DNA=record.seq

# Printing out positions of GRINS regions in the initial sequence, and creating an annotated GenBank file
if len(region_starts)!=len(region_ends):
    GRINS_file.write("Problem: non-equal number of GRINS start and end locations")
    print "Problem: non-equal number of GRINS start and end locations"
elif len(region_starts)==0 and len(region_ends)==0 and len(possible_region_starts)==0 and len(possible_region_ends)==0:
    print "No GRINS found in NIDS"
    GRINS_file.write("No GRINS found in NIDS")
else:
    GRINS_file.write("With cutoff values for GRINS start above "+str(local_start_cutoff)+" for local identity and above "+str(next10_start_cutoff)+" for the next 10 windows\n")
    GRINS_file.write("and cutoff values for GRINS end below "+str(local_end_cutoff)+" for local identity and above "+str(previous10_end_cutoff)+" for the previous 10 windows,\n")
    GRINS_file.write("found following GRINS in NIDS:\n")
with open("NIDS_GRINS.gb","w") as annotated_file:
    record_annotated=record
    for i in range(0,len(region_starts)):
        GRINS_file.write(str(region_starts[i]*30)+"-"+str(region_ends[i]*30)+"\n")
        feature=SeqFeature(FeatureLocation(start=region_starts[i]*30, end=region_ends[i]*30), type='GRINS')
        record_annotated.features.append(feature)
    for i in range(0,len(possible_region_starts)):
        GRINS_file.write(str(possible_region_starts[i]*30)+"-"+str(possible_region_ends[i]*30)+" (erased/possible GRINS)\n")
        feature=SeqFeature(FeatureLocation(start=possible_region_starts[i]*30, end=possible_region_ends[i]*30), type='GRINS?')
        record_annotated.features.append(feature)
    SeqIO.write(record_annotated, annotated_file, 'genbank')


# Constructing lists of skew values for the entire sequence
GCskew_array=[]
TAskew_array=[]
ATskew_array=[]
GCcontent_array=[]
A_number_array=[]
T_number_array=[]
G_number_array=[]
C_number_array=[]
for i in range(0,len(input_DNA)-149,30):
    fragment=input_DNA[i:i+150]
    G_number,C_number,GCskew,GCcontent=GC_skew(fragment)
    A_number,T_number,TAskew=TA_skew(fragment)
    A_number,T_number,ATskew=AT_skew(fragment)
    GCskew_array.append(GCskew)
    TAskew_array.append(TAskew)
    ATskew_array.append(ATskew)
    GCcontent_array.append(GCcontent)
    A_number_array.append(float(A_number)/150.0)
    T_number_array.append(float(T_number)/150.0)
    G_number_array.append(float(G_number)/150.0)
    C_number_array.append(float(C_number)/150.0)

a1=range(0,len(GCskew_array))
a2=range(0,len(similarity_max_array))

### The following section is only executed if there are GRINSs in the sequence: ###

if len(region_starts)!=0 and len(region_ends)!=0:
    
    # Constructing lists of skew values for GRINS and non-GRINS sequences, if they exist
    A_number_array_GRINS=[]
    G_number_array_GRINS=[]
    C_number_array_GRINS=[]
    T_number_array_GRINS=[]
    GCskew_array_GRINS=[]
    TAskew_array_GRINS=[]
    ATskew_array_GRINS=[]
    GCskew_array_GRINS_absolute=[]
    TAskew_array_GRINS_absolute=[]
    ATskew_array_GRINS_absolute=[]
    GCcontent_array_GRINS=[]
    A_number_array_rest=[]
    G_number_array_rest=[]
    C_number_array_rest=[]
    T_number_array_rest=[]
    GCskew_array_rest=[]
    TAskew_array_rest=[]
    ATskew_array_rest=[]
    GCskew_array_rest_absolute=[]
    TAskew_array_rest_absolute=[]
    ATskew_array_rest_absolute=[]
    GCcontent_array_rest=[]
    curr_end=0
    for i in range(0,len(region_starts)):
        for j in range(curr_end,region_starts[i]):
            A_number_array_rest.append(A_number_array[j])
            G_number_array_rest.append(G_number_array[j])
            C_number_array_rest.append(C_number_array[j])
            T_number_array_rest.append(T_number_array[j]) 
            GCskew_array_rest.append(GCskew_array[j])  
            TAskew_array_rest.append(TAskew_array[j])
            ATskew_array_rest.append(ATskew_array[j])
            GCskew_array_rest_absolute.append(abs(GCskew_array[j]))  
            TAskew_array_rest_absolute.append(abs(TAskew_array[j]))
            ATskew_array_rest_absolute.append(abs(ATskew_array[j]))
            GCcontent_array_rest.append(GCcontent_array[j]) 
        for k in range(region_starts[i],region_ends[i]):
            A_number_array_GRINS.append(A_number_array[k])
            G_number_array_GRINS.append(G_number_array[k])
            C_number_array_GRINS.append(C_number_array[k])
            T_number_array_GRINS.append(T_number_array[k])
            GCskew_array_GRINS.append(GCskew_array[k])  
            TAskew_array_GRINS.append(TAskew_array[k])
            ATskew_array_GRINS.append(ATskew_array[k])
            GCskew_array_GRINS_absolute.append(abs(GCskew_array[k]))  
            TAskew_array_GRINS_absolute.append(abs(TAskew_array[k]))
            ATskew_array_GRINS_absolute.append(abs(ATskew_array[k]))
            GCcontent_array_GRINS.append(GCcontent_array[k]) 
        curr_end=region_ends[i]        
    for j in range(curr_end,len(A_number_array)):
        A_number_array_rest.append(A_number_array[j])
        G_number_array_rest.append(G_number_array[j])
        C_number_array_rest.append(C_number_array[j])
        T_number_array_rest.append(T_number_array[j])
        GCskew_array_rest.append(GCskew_array[j])  
        TAskew_array_rest.append(TAskew_array[j])
        ATskew_array_rest.append(ATskew_array[j])
        GCcontent_array_rest.append(GCcontent_array[j]) 

    # Transforming lists into arrays
    G_number_array_GRINS2=np.array(G_number_array_GRINS)
    C_number_array_GRINS2=np.array(C_number_array_GRINS)
    A_number_array_GRINS2=np.array(A_number_array_GRINS)
    T_number_array_GRINS2=np.array(T_number_array_GRINS)
    G_number_array_rest2=np.array(G_number_array_rest)
    C_number_array_rest2=np.array(C_number_array_rest)
    A_number_array_rest2=np.array(A_number_array_rest)
    T_number_array_rest2=np.array(T_number_array_rest)
    GCskew_array_GRINS2=np.array(GCskew_array_GRINS)
    TAskew_array_GRINS2=np.array(TAskew_array_GRINS)
    ATskew_array_GRINS2=np.array(ATskew_array_GRINS)
    GCskew_array_rest2=np.array(GCskew_array_rest)
    TAskew_array_rest2=np.array(TAskew_array_rest)
    ATskew_array_rest2=np.array(ATskew_array_rest)


### Plotting graphs ###


## Plotting graphs for skews in GRINS and non-GRINS regions and max identity scores: ##
# GC skews are in blue
# AT skews are in red
# TA skews are in green
# Baseline at 0 is in gray
# Max identity scores are in black
# GRINS regions are shaded gray

# GC and TA skews; similarity
plt.figure(1)
fig,ax1=plt.subplots(figsize=(10, 4))
ax1.set_ylim(-1,1)
ax1.set_xlim(0,len(GCskew_array))
ax1.grid(False)
plt.plot(a1,GCskew_array,color="blue",linewidth=0.5)
plt.plot(a1,TAskew_array,color="green",linewidth=0.5)
for i in range(0,len(region_ends)):
    rect = patches.Rectangle((region_starts[i],1),region_ends[i]-region_starts[i],-2,linewidth=1,edgecolor='gray',facecolor='gray', alpha=0.3)
    ax1.add_patch(rect)
for i in range(0,len(possible_region_ends)):
    rect = patches.Rectangle((possible_region_starts[i],1),possible_region_ends[i]-possible_region_starts[i],-2,linewidth=1,edgecolor='lightgray',facecolor='lightgray', alpha=0.3)
    ax1.add_patch(rect)
plt.legend(["GC skew","TA skew"],loc=3)
ax2=ax1.twinx()
ax2.set_ylim(0,100)
ax2.set_xlim(0,len(similarity_max_array))
ax2.grid(False)
plt.plot(a2,similarity_max_array,color="black",linewidth=0.5)
plt.legend(["Maximum identity %"],loc=4)
plt.savefig("NIDS_GCskew_TAskew_identity.png")

# GC and AT skews; similarity
plt.figure(2)
fig,ax1=plt.subplots(figsize=(10, 4))
ax1.set_ylim(-1,1)
ax1.set_xlim(0,len(GCskew_array))
ax1.grid(False)
plt.plot(a1,GCskew_array,color="blue",linewidth=0.5)
plt.plot(a1,ATskew_array,color="red",linewidth=0.5)
for i in range(0,len(region_ends)):
    rect = patches.Rectangle((region_starts[i],1),region_ends[i]-region_starts[i],-2,linewidth=1,edgecolor='gray',facecolor='gray', alpha=0.3)
    ax1.add_patch(rect)
for i in range(0,len(possible_region_ends)):
    rect = patches.Rectangle((possible_region_starts[i],1),possible_region_ends[i]-possible_region_starts[i],-2,linewidth=1,edgecolor='lightgray',facecolor='lightgray', alpha=0.3)
    ax1.add_patch(rect)
plt.legend(["GC skew","AT skew"],loc=3)
ax2=ax1.twinx()
ax2.set_ylim(0,100)
ax2.set_xlim(0,len(similarity_max_array))
ax2.grid(False)
plt.plot(a2,similarity_max_array,color="black",linewidth=0.5)
plt.legend(["Maximum identity %"],loc=4)
plt.savefig("NIDS_GCskew_ATskew_identity.png")


if len(region_starts)!=0 and len(region_ends)!=0:

    ## Plotting correlation graphs for all nucleotide pairs in GRINS and non-GRINS regions, and printing the correlation coefficients, only if they exist ##
    os.mkdir("./Correlation_nt")
    skew_corr_file=open("NIDS_composition_and_skew_values.txt",'w')

    # Correlation of G and C nucleotides
    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(G_number_array_GRINS, C_number_array_GRINS)
    skew_corr_file.write("Correlation between G and C in GRINS: r="+str(np.around(r_value,2))+", p="+str(np.around(p_value,2))+"\n")
    plt.figure(3)
    fig,ax1=plt.subplots(figsize=(4, 4))
    ax1.set_ylim(0,0.7)
    ax1.set_xlim(0,0.7)
    ax1.text(0.4,0.55,"r="+str(np.round(r_value,2))+p_value_def(p_value), fontdict=font)
    sns.regplot(x=G_number_array_GRINS2,y=C_number_array_GRINS2, fit_reg=True)
    plt.savefig("./Correlation_nt/Correlation_G_C_NIDS_GRINS.png")

    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(G_number_array_rest, C_number_array_rest)
    skew_corr_file.write("Correlation between G and C in non-GRINS regions: r="+str(np.around(r_value,2))+", p="+str(np.around(p_value,2))+"\n")
    plt.figure(4)
    fig,ax1=plt.subplots(figsize=(4, 4))
    ax1.set_ylim(0,0.7)
    ax1.set_xlim(0,0.7)
    ax1.text(0.4,0.55,"r="+str(np.round(r_value,2))+p_value_def(p_value), fontdict=font)
    sns.regplot(x=G_number_array_rest2,y=C_number_array_rest2, fit_reg=True)
    plt.savefig("./Correlation_nt/Correlation_G_C_NIDS_rest.png")

    # Correlation of G and A nucleotides
    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(G_number_array_GRINS, A_number_array_GRINS)
    skew_corr_file.write("Correlation between G and A in GRINS: r="+str(np.around(r_value,2))+", p="+str(np.around(p_value,2))+"\n")
    plt.figure(5)
    fig,ax1=plt.subplots(figsize=(4, 4))
    ax1.set_ylim(0,0.7)
    ax1.set_xlim(0,0.7)
    ax1.text(0.4,0.55,"r="+str(np.round(r_value,2))+p_value_def(p_value), fontdict=font)
    sns.regplot(x=G_number_array_GRINS2,y=A_number_array_GRINS2, fit_reg=True)
    plt.savefig("./Correlation_nt/Correlation_G_A_NIDS_GRINS.png")

    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(G_number_array_rest, A_number_array_rest)
    skew_corr_file.write("Correlation between G and A in non-GRINS regions: r="+str(np.around(r_value,2))+", p="+str(np.around(p_value,2))+"\n")
    plt.figure(6)
    fig,ax1=plt.subplots(figsize=(4, 4))
    ax1.set_ylim(0,0.7)
    ax1.set_xlim(0,0.7)
    ax1.text(0.4,0.55,"r="+str(np.round(r_value,2))+p_value_def(p_value), fontdict=font)
    sns.regplot(x=G_number_array_rest2,y=A_number_array_rest2, fit_reg=True)
    plt.savefig("./Correlation_nt/Correlation_G_A_NIDS_rest.png")

    # Correlation of G and T nucleotides
    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(G_number_array_GRINS, T_number_array_GRINS)
    skew_corr_file.write("Correlation between G and T in GRINS: r="+str(np.around(r_value,2))+", p="+str(np.around(p_value,2))+"\n")
    plt.figure(7)
    fig,ax1=plt.subplots(figsize=(4, 4))
    ax1.set_ylim(0,0.7)
    ax1.set_xlim(0,0.7)
    ax1.text(0.4,0.55,"r="+str(np.round(r_value,2))+p_value_def(p_value), fontdict=font)
    sns.regplot(x=G_number_array_GRINS2,y=T_number_array_GRINS2, fit_reg=True)
    plt.savefig("./Correlation_nt/Correlation_G_T_NIDS_GRINS.png")

    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(G_number_array_rest, T_number_array_rest)
    skew_corr_file.write("Correlation between G and T in non-GRINS regions: r="+str(np.around(r_value,2))+", p="+str(np.around(p_value,2))+"\n")
    plt.figure(8)
    fig,ax1=plt.subplots(figsize=(4, 4))
    ax1.set_ylim(0,0.7)
    ax1.set_xlim(0,0.7)
    ax1.text(0.4,0.55,"r="+str(np.round(r_value,2))+p_value_def(p_value), fontdict=font)
    sns.regplot(x=G_number_array_rest2,y=T_number_array_rest2, fit_reg=True)
    plt.savefig("./Correlation_nt/Correlation_G_T_NIDS_rest.png")

    # Correlation of A and C nucleotides
    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(A_number_array_GRINS, C_number_array_GRINS)
    skew_corr_file.write("Correlation between A and C in GRINS: r="+str(np.around(r_value,2))+", p="+str(np.around(p_value,2))+"\n")
    plt.figure(9)
    fig,ax1=plt.subplots(figsize=(4, 4))
    ax1.set_ylim(0,0.7)
    ax1.set_xlim(0,0.7)
    ax1.text(0.4,0.55,"r="+str(np.round(r_value,2))+p_value_def(p_value), fontdict=font)
    sns.regplot(x=A_number_array_GRINS2,y=C_number_array_GRINS2, fit_reg=True)
    plt.savefig("./Correlation_nt/Correlation_A_C_NIDS_GRINS.png")

    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(A_number_array_rest, C_number_array_rest)
    skew_corr_file.write("Correlation between A and C in non-GRINS regions: r="+str(np.around(r_value,2))+", p="+str(np.around(p_value,2))+"\n")
    plt.figure(10)
    fig,ax1=plt.subplots(figsize=(4, 4))
    ax1.set_ylim(0,0.7)
    ax1.set_xlim(0,0.7)
    ax1.text(0.4,0.55,"r="+str(np.round(r_value,2))+p_value_def(p_value), fontdict=font)
    sns.regplot(x=A_number_array_rest2,y=C_number_array_rest2, fit_reg=True)
    plt.savefig("./Correlation_nt/Correlation_A_C_NIDS_rest.png")

    # Correlation of C and T nucleotides
    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(C_number_array_GRINS, T_number_array_GRINS)
    skew_corr_file.write("Correlation between C and T in GRINS: r="+str(np.around(r_value,2))+", p="+str(np.around(p_value,2))+"\n")
    plt.figure(11)
    fig,ax1=plt.subplots(figsize=(4, 4))
    ax1.set_ylim(0,0.7)
    ax1.set_xlim(0,0.7)
    ax1.text(0.4,0.55,"r="+str(np.round(r_value,2))+p_value_def(p_value), fontdict=font)
    sns.regplot(x=C_number_array_GRINS2,y=T_number_array_GRINS2, fit_reg=True)
    plt.savefig("./Correlation_nt/Correlation_C_T_NIDS_GRINS.png")

    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(C_number_array_rest, T_number_array_rest)
    skew_corr_file.write("Correlation between C and T in non-GRINS regions: r="+str(np.around(r_value,2))+", p="+str(np.around(p_value,2))+"\n")
    plt.figure(12)
    fig,ax1=plt.subplots(figsize=(4, 4))
    ax1.set_ylim(0,0.7)
    ax1.set_xlim(0,0.7)
    ax1.text(0.4,0.55,"r="+str(np.round(r_value,2))+p_value_def(p_value), fontdict=font)
    sns.regplot(x=C_number_array_rest2,y=T_number_array_rest2, fit_reg=True)
    plt.savefig("./Correlation_nt/Correlation_C_T_NIDS_rest.png")

    # Correlation of A and T nucleotides
    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(A_number_array_GRINS, T_number_array_GRINS)
    skew_corr_file.write("Correlation between A and T in GRINS: r="+str(np.around(r_value,2))+", p="+str(np.around(p_value,2))+"\n")
    plt.figure(13)
    fig,ax1=plt.subplots(figsize=(4, 4))
    ax1.set_ylim(0,0.7)
    ax1.set_xlim(0,0.7)
    ax1.text(0.4,0.55,"r="+str(np.round(r_value,2))+p_value_def(p_value), fontdict=font)
    sns.regplot(x=A_number_array_GRINS2,y=T_number_array_GRINS2, fit_reg=True)
    plt.savefig("./Correlation_nt/Correlation_A_T_NIDS_GRINS.png")

    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(A_number_array_rest, T_number_array_rest)
    skew_corr_file.write("Correlation between A and T in non-GRINS regions: r="+str(np.around(r_value,2))+", p="+str(np.around(p_value,2))+"\n")
    plt.figure(14)
    fig,ax1=plt.subplots(figsize=(4, 4))
    ax1.set_ylim(0,0.7)
    ax1.set_xlim(0,0.7)
    ax1.text(0.4,0.55,"r="+str(np.round(r_value,2))+p_value_def(p_value), fontdict=font)
    sns.regplot(x=A_number_array_rest2,y=T_number_array_rest2, fit_reg=True)
    plt.savefig("./Correlation_nt/Correlation_A_T_NIDS_rest.png")


    ## Plotting correlation graphs for SKEWS in GRINS and non-GRINS regions, and printing the correlation coefficients ##
    os.mkdir("./Correlation_skews")

    # Correlation of GC and TA skews
    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(GCskew_array_GRINS, TAskew_array_GRINS)
    skew_corr_file.write("Correlation between GC and TA skews in GRINS: r="+str(np.around(r_value,2))+", p="+str(np.around(p_value,2))+"\n")
    plt.figure(15)
    fig,ax1=plt.subplots(figsize=(4, 4))
    ax1.set_ylim(-1.1,1.1)
    ax1.set_xlim(-1.1,1.1)
    sns.regplot(x=GCskew_array_GRINS2,y=TAskew_array_GRINS2, fit_reg=True)
    ax1.text(-1,0.6,"r="+str(np.round(r_value,2))+p_value_def(p_value), fontdict=font)
    plt.savefig("./Correlation_skews/Correlation_GC_TA_skews_NIDS_GRINS.png")

    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(GCskew_array_rest, TAskew_array_rest)
    skew_corr_file.write("Correlation between GC and TA skews in non-GRINS regions: r="+str(np.around(r_value,2))+", p="+str(np.around(p_value,2))+"\n")
    plt.figure(16)
    fig,ax1=plt.subplots(figsize=(4, 4))
    ax1.set_ylim(-1.1,1.1)
    ax1.set_xlim(-1.1,1.1)
    sns.regplot(x=GCskew_array_rest2,y=TAskew_array_rest2, fit_reg=True)
    ax1.text(-1,0.6,"r="+str(np.round(r_value,2))+p_value_def(p_value), fontdict=font)
    plt.savefig("./Correlation_skews/Correlation_GC_TA_skews_NIDS_rest.png")

    # Correlation of GC and AT skews
    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(GCskew_array_GRINS, ATskew_array_GRINS)
    skew_corr_file.write("Correlation between GC and AT skews in GRINS: r="+str(np.around(r_value,2))+", p="+str(np.around(p_value,2))+"\n")   
    plt.figure(17)
    fig,ax1=plt.subplots(figsize=(4, 4))
    ax1.set_ylim(-1.1,1.1)
    ax1.set_xlim(-1.1,1.1)
    sns.regplot(x=GCskew_array_GRINS2,y=ATskew_array_GRINS2, fit_reg=True)
    ax1.text(-1,0.6,"r="+str(np.round(r_value,2))+p_value_def(p_value), fontdict=font)
    plt.savefig("./Correlation_skews/Correlation_GC_AT_skews_NIDS_GRINS.png")

    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(GCskew_array_rest, ATskew_array_rest)
    skew_corr_file.write("Correlation between GC and AT skews in non-GRINS regions: r="+str(np.around(r_value,2))+", p="+str(np.around(p_value,2))+"\n"+"\n")   
    plt.figure(18)
    fig,ax1=plt.subplots(figsize=(4, 4))
    ax1.set_ylim(-1.1,1.1)
    ax1.set_xlim(-1.1,1.1)
    sns.regplot(x=GCskew_array_rest2,y=ATskew_array_rest2, fit_reg=True)
    ax1.text(-1,0.6,"r="+str(np.round(r_value,2))+p_value_def(p_value), fontdict=font)
    plt.savefig("./Correlation_skews/Correlation_GC_AT_skews_NIDS_rest.png")
    

    ## Saving mean, standard deviation, min and max of absolute skew values into file
    skew_corr_file.write("Mean absolute GC skew in GRINSs: "+str(np.mean(GCskew_array_GRINS_absolute))+"\n")
    skew_corr_file.write("Deviation of mean absolute GC skew in GRINSs: "+str(np.std(GCskew_array_GRINS_absolute))+"\n")
    skew_corr_file.write("Maximum absolute GC skew in GRINSs: "+str(np.max(GCskew_array_GRINS_absolute))+"\n")
    skew_corr_file.write("Minimum absolute GC skew in GRINSs: "+str(np.min(GCskew_array_GRINS_absolute))+"\n")

    skew_corr_file.write("Mean absolute TA skew in GRINSs: "+str(np.mean(TAskew_array_GRINS_absolute))+"\n")
    skew_corr_file.write("Deviation of mean absolute TA skew in GRINSs: "+str(np.std(TAskew_array_GRINS_absolute))+"\n")
    skew_corr_file.write("Maximum absolute TA skew in GRINSs: "+str(np.max(TAskew_array_GRINS_absolute))+"\n")
    skew_corr_file.write("Minimum absolute TA skew in GRINSs: "+str(np.min(TAskew_array_GRINS_absolute))+"\n")

    skew_corr_file.write("Mean absolute GC skew in non-GRINS regions: "+str(np.mean(GCskew_array_rest_absolute))+"\n")
    skew_corr_file.write("Deviation of mean absolute GC skew in non-GRINS regions: "+str(np.std(GCskew_array_rest_absolute))+"\n")
    skew_corr_file.write("Maximum absolute GC skew in non-GRINS regions: "+str(np.max(GCskew_array_rest_absolute))+"\n")
    skew_corr_file.write("Minimum absolute GC skew in non-GRINS regions: "+str(np.min(GCskew_array_rest_absolute))+"\n")

    skew_corr_file.write("Mean absolute TA skew in non-GRINS regions:: "+str(np.mean(TAskew_array_rest_absolute))+"\n")
    skew_corr_file.write("Deviation of mean absolute TA skew in non-GRINS regions: "+str(np.std(TAskew_array_rest_absolute))+"\n")
    skew_corr_file.write("Maximum absolute TA skew in non-GRINS regions: "+str(np.max(TAskew_array_rest_absolute))+"\n")
    skew_corr_file.write("Minimum absolute ta skew in non-GRINS regions: "+str(np.min(TAskew_array_rest_absolute))+"\n")

## Plotting nucleotide content graphs ##
os.mkdir("./Content_nt")

# G and C content
plt.figure(19)
fig,ax1=plt.subplots(figsize=(10, 4))
ax1.set_ylim(0,1)
ax1.set_xlim(0,len(A_number_array))
ax1.grid(False)
plt.plot(a1,G_number_array,color="orange",linewidth=0.5)
plt.plot(a1,C_number_array,color="blue",linewidth=0.5)
for i in range(0,len(region_ends)):
    rect = patches.Rectangle((region_starts[i],1),region_ends[i]-region_starts[i],-2,linewidth=1,edgecolor='gray',facecolor='gray', alpha=0.3)
    ax1.add_patch(rect)
for i in range(0,len(possible_region_ends)):
    rect = patches.Rectangle((possible_region_starts[i],1),possible_region_ends[i]-possible_region_starts[i],-2,linewidth=1,edgecolor='lightgray',facecolor='lightgray', alpha=0.3)
    ax1.add_patch(rect)
plt.legend(["G","C"],loc=1)
plt.savefig("./Content_nt/NIDS_nucleotide_G_C_content.png")

# A and T content
plt.figure(20)
fig,ax1=plt.subplots(figsize=(10, 4))
ax1.set_ylim(0,1)
ax1.set_xlim(0,len(A_number_array))
ax1.grid(False)
plt.plot(a1,A_number_array,color="red",linewidth=0.5)
plt.plot(a1,T_number_array,color="green",linewidth=0.5)
for i in range(0,len(region_ends)):
    rect = patches.Rectangle((region_starts[i],1),region_ends[i]-region_starts[i],-2,linewidth=1,edgecolor='gray',facecolor='gray', alpha=0.3)
    ax1.add_patch(rect)
for i in range(0,len(possible_region_ends)):
    rect = patches.Rectangle((possible_region_starts[i],1),possible_region_ends[i]-possible_region_starts[i],-2,linewidth=1,edgecolor='lightgray',facecolor='lightgray', alpha=0.3)
    ax1.add_patch(rect)
plt.legend(["A","T"],loc=1)
plt.savefig("./Content_nt/NIDS_nucleotide_A_T_content.png")

# GC content
plt.figure(21)
fig,ax1=plt.subplots(figsize=(10, 4))
ax1.set_ylim(0,1)
ax1.set_xlim(0,len(GCskew_array))
ax1.grid(False)
plt.plot(a1,GCcontent_array,color="black",linewidth=0.5)
for i in range(0,len(region_ends)):
    rect = patches.Rectangle((region_starts[i],1),region_ends[i]-region_starts[i],-2,linewidth=1,edgecolor='gray',facecolor='gray', alpha=0.3)
    ax1.add_patch(rect)
for i in range(0,len(possible_region_ends)):
    rect = patches.Rectangle((possible_region_starts[i],1),possible_region_ends[i]-possible_region_starts[i],-2,linewidth=1,edgecolor='lightgray',facecolor='lightgray', alpha=0.3)
    ax1.add_patch(rect)
plt.legend(["GC content"],loc=4)
plt.savefig("./Content_nt/NIDS_GC_content.png")


plt.figure(22)
fig,ax1=plt.subplots(figsize=(10, 4))
ax1.set_ylim(-1,1)
ax1.set_xlim(300,550)
ax1.grid(False)
plt.plot(a1,GCskew_array,color="blue",linewidth=0.5)
plt.plot(a1,TAskew_array,color="green",linewidth=0.5)
for i in range(0,len(region_ends)):
    rect = patches.Rectangle((region_starts[i],1),region_ends[i]-region_starts[i],-2,linewidth=1,edgecolor='gray',facecolor='gray', alpha=0.3)
    ax1.add_patch(rect)
for i in range(0,len(possible_region_ends)):
    rect = patches.Rectangle((possible_region_starts[i],1),possible_region_ends[i]-possible_region_starts[i],-2,linewidth=1,edgecolor='lightgray',facecolor='lightgray', alpha=0.3)
    ax1.add_patch(rect)
plt.legend(["GC skew","TA skew"],loc=3)
ax2=ax1.twinx()
ax2.set_ylim(0,100)
ax2.set_xlim(300,550)
ax2.grid(False)
plt.plot(a2,similarity_max_array,color="black",linewidth=0.5)
plt.legend(["Maximum identity %"],loc=4)
plt.savefig("NIDS_GCskew_TAskew_identity_segment.png")


plt.figure(23)
fig,ax1=plt.subplots(figsize=(10, 4))
ax1.set_ylim(-1,1)
ax1.set_xlim(0,len(GCskew_array))
ax1.grid(False)
plt.plot(a1,GCskew_array,color="blue",linewidth=0.5)
plt.plot(a1,TAskew_array,color="green",linewidth=0.5)
plt.legend(["GC skew","TA skew"],loc=3)
plt.savefig("NIDS_GCskew_TAskew.png")


