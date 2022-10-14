import glob, os, sys, argparse
import pandas as pd
#import seaborn as sns
#import matplotlib.pyplot as plt
#import matplotlib
#matplotlib.rcParams['pdf.fonttype'] = 42
#matplotlib.rcParams['ps.fonttype'] = 42

from collections import defaultdict, namedtuple

## DEFINE SOME VARIABLES WE NEED: kreport files, extended taxon search terms, theoretical filenames, etc

# Some taxa need modified search terms. This dict should work for all datasets.
taxa_searchterms = {'Limosilactobacillus fermentum':'Lactobacillus fermentum|Limosilactobacillus fermentum',
                       "Veillonella rogosae": "Veillonella",
                      "Prevotella corporis": "Prevotella",
                      'Bacillus subtilis': 'Bacillus subtilis|Bacillus spizizenii',
                      'Luteovulum sphaeroides':'Rhodobacter sphaeroides|Luteovulum sphaeroides|Cereibacter sphaeroides',
                      'Phocaeicola vulgatus': 'Bacteroides vulgatus|Phocaeicola vulgatus',
                      'Luteovulum':'Rhodobacter|Luteovulum|Cereibacter',
                      'Phocaeicola':'Bacteroides|Phocaeicola',
                      'Limosilactobacillus': 'Lactobacillus|Limosilactobacillus'
                       }

all_datasets = {"ATCC_MSA1003": ["HiFi_ATCC_MSA1003_full", "HiFi_ATCC_MSA1003_short", "Illumina_ATCC_MSA1003"],
              "Zymo_D6300": ["Illumina_Zymo_D6300","ONT_Q20_Zymo_D6300_full","ONT_Q20_Zymo_D6300_short","ONT_R10_Zymo_D6300"],
              "HiFi_Zymo_D6331": ["HiFi_Zymo_D6331"]}

tools= ["Bracken", "BugSeq-V2", "Centrifuge-h22",
        "Centrifuge-h500", "Kraken2", "MEGAN-LR-Nuc-HiFi",
        "MEGAN-LR-Nuc-ONT", "MEGAN-LR-Prot", "MMseqs2",
        "Metamaps", "Metaphlan3", "mOTUs"]

# filenames for all theoretical distributions
species_filenames = {"Zymo_D6300": "Zymo_D6300.theoretical-distrib.species.csv",
                     "HiFi_Zymo_D6331": "HiFi_Zymo_D6331.theoretical-distrib.species.csv",
                     "ATCC_MSA1003": "ATCC_MSA1003.theoretical-distrib.species.csv"}

genus_filenames = {"Zymo_D6300": "Zymo_D6300.theoretical-distrib.genus.csv",
                     "HiFi_Zymo_D6331": "HiFi_Zymo_D6331.theoretical-distrib.genus.csv",
                     "ATCC_MSA1003": "ATCC_MSA1003.theoretical-distrib.genus.csv"}

theoretical_distributions_dict = {"species": species_filenames, "genus": genus_filenames}

# filenames for all sourmash kreports
df_thresh_folder = "sourmash-kreports-2022-10-03.genbank"
zero_thresh_folder = "../output.lr-gather/gather"

kreport_dict = defaultdict(dict)

for kreport_folder in [df_thresh_folder]:#, zero_thresh_folder]:
    for ksize in [31,51]:
        for ds, ds_list in all_datasets.items():
            for d in ds_list:
                fname = f"Sourmash-k{ksize}"
                kreportF = f"{kreport_folder}/{d}.dna-k{ksize}-sc1000.gather.genbank.kreport.txt"
                kreport_dict[d].update({fname: kreportF})
for kreport_folder in [zero_thresh_folder]:
    for ksize in [31,51]:
        for ds, ds_list in all_datasets.items():
            for d in ds_list:
                fname = f"Sourmash-k{ksize}_zeroThresh"
                kreportF = f"{kreport_folder}/{d}.dna-k{ksize}-sc1000.gather.genbank.kreport.txt"
                kreport_dict[d].update({fname: kreportF})

# now add non-sourmash kreports
full_kreport_files = glob.glob("All-kreports/*/inputs/*.kreport.txt")
for f in full_kreport_files:
    if "D6331" in f:
        rd = "HiFi_Zymo_D6331"
    elif "ATCC" in f:
        if "Illumina" in f:
            rd = "Illumina_ATCC_MSA1003"
        else:
            rd = "HiFi_ATCC_MSA1003_full"
    else:
        if "R10" in f:
            rd = "ONT_R10_Zymo_D6300"
        elif "Illumina" in f:
            rd = "Illumina_Zymo_D6300"
        else:
            rd = 'ONT_Q20_Zymo_D6300_full'
    tool = os.path.basename(f).split('.')[0]
    if "Sourmash" in tool: # no need to add these twice
        continue
    kreport_dict[rd].update({tool:f})

# ADD kreports for SIMULATED SHORT  READS
sim_kreport_files = glob.glob("All-kreports/simread-kreports/*.kreport.txt")
#print(len(sim_kreport_files))
for f in sim_kreport_files:
    if "ATCC" in f:
        rd = "HiFi_ATCC_MSA1003_short"
    elif "D6300" in f:
        rd = 'ONT_Q20_Zymo_D6300_short'
    tool = os.path.basename(f).split('.')[0]
    if "Sourmash" in tool: # no need to add these twice
        continue
    kreport_dict[rd].update({tool:f})

# define taxresult namedtuple
TaxResult = namedtuple('TaxResult', "taxon,rank,mock_community,read_dataset,tool,d1perc,dp1perc,dp01perc,dany,cumulative_count")
    # set minimum detection threshold for this read dataset
    #detect_dict = {'0.001-perc': int(0.00001*total_reads),
    #           '0.1-perc': int(0.001*total_reads),
    #           '1-perc': int(0.01*total_reads)}

    # change this to generate results for other filtering threshold
   # dlevel = '0.001-perc'

    #for k, v in detect_dict.items():
    #    print(k, '=', v)

def assess_kreport(dataset_name, read_dataset, theoretical_dict, kreport_fname, kreport_df, namediff_dict, detection_thresholds = [0.01, 0.001, 0.0001, 0]):
    """
    Read in a kreport from a read dataset; parse and store results
    """
    taxresults = []
    rank_labels = {"species": 'S', "genus": 'G'}

    # First, search for True Positives and False Negatives
    for rank in rank_labels.keys():
        rk = rank_labels[rank]
        theoretical_info = theoretical_dict[rank]
        theoretical_taxa = theoretical_info['Taxon'].tolist()
        searchnames = []
        # loop through expected taxa. First grab theoretical/expected reads/percent, use to assess tool accuracy
        for t in theoretical_taxa:
            reads_column = f"{read_dataset}_theoretical_reads"
            theoretical_reads = theoretical_info[theoretical_info['Taxon'] == t][reads_column].values[0]
            # Now, look in the kreport for this taxon, using the extended taxon searchname, if there is one
            searchname = namediff_dict.get(t, t)
            searchnames.append(searchname)
            if "|" in searchname:
                tempD = kreport_df[(kreport_df['name'].str.contains(searchname, regex=True))]
            else:
                tempD = kreport_df[(kreport_df['name'].str.contains(t))]
            # find the cumulative count for the taxon
            # for species, this will have the effect of including all strains below, which are correct
            this_count = tempD.loc[kreport_df['rank'] == rk, 'cumulative_count'].sum()
            detection_info = []
            # for sourmash, check if % exceeds minimum % detection threshold
            if "Sourmash" in kreport_fname:
                theoretical_percent = theoretical_info[theoretical_info['Taxon'] == t]["Theoretical Distribution"].values[0]
                this_perc = tempD.loc[kreport_df['rank'] == rk, 'proportion'].sum() # same idea as for counts
                # now, check results at the different detection thresholds.
                for thresh in detection_thresholds:
                    res ="FN"
                    perc_needed = theoretical_percent * thresh
                    if this_perc >= perc_needed and this_perc >0:
                        res="TP"
                    detection_info.append(res)
            # for all other tools, check if read count exceeds minimum read threshold
            else:
                for thresh in detection_thresholds:
                    res ="FN"
                    reads_needed = theoretical_reads * thresh
                    #print(this_count, reads_needed)
                    #import pdb;pdb.set_trace()
                    if this_count >= reads_needed and this_count !=0:
                        res="TP"
                    detection_info.append(res)
            taxreslist = [t, rank, dataset_name, read_dataset, kreport_fname] + detection_info + [this_count]
            taxres = TaxResult(*taxreslist)
            taxresults.append(taxres)

        # Now, look for False Positives: eliminate all of expected taxa, then loop through remaining and record
        fp_df= kreport_df[(~kreport_df['name'].str.contains("|".join(searchnames), regex=True))]
        rank_fp_df =  fp_df.loc[fp_df['rank'] == rk]
        for taxon in rank_fp_df['name'].to_list():
            this_taxon_count = rank_fp_df[rank_fp_df["name"] == taxon]['cumulative_count'].sum()
            taxres = TaxResult(taxon, rank, dataset_name, read_dataset, kreport_fname, "FP","FP","FP","FP", this_taxon_count)
            taxresults.append(taxres)
    return taxresults


def main(args):
    all_taxres = []

    for mock_community, read_datasets in all_datasets.items():
        for rd in read_datasets:
            # get theoretical distribution dataframes for this read dataset
            theoretical_species = pd.read_csv(theoretical_distributions_dict["species"][mock_community])
            theoretical_genera = pd.read_csv(theoretical_distributions_dict["genus"][mock_community])
            theoretical = {"species": theoretical_species, "genus": theoretical_genera}
            # get results for each kreport
            for kreport_fname, kreportF in kreport_dict[rd].items(): # for tool, kreportfile ...
                # read kreport
                kreport_df = pd.read_csv(kreportF, sep = '\t', header=None,
                         names = ['proportion', 'cumulative_count', 'level_count', 'rank', 'taxid', 'name'])
                # assess vs theoretical taxon distributions
                taxres = assess_kreport(mock_community, rd, theoretical, kreport_fname, kreport_df, taxa_searchterms)
                all_taxres +=taxres

    # now convert all the results to a dataframe
    all_results = pd.DataFrame.from_records(all_taxres, columns=TaxResult._fields)
    all_results.head()
    # write to csv
    all_results.to_csv("kreport-results.csv.gz", index=False)

def cmdline(sys_args):
    "Command line entry point w/argparse action."
    p = argparse.ArgumentParser()
#    p.add_argument("--file-dir", default="notebooks")
    args = p.parse_args()
    return main(args)

if __name__ == '__main__':
    returncode = cmdline(sys.argv[1:])
    sys.exit(returncode)


# some tests to see how it's going!

def make_theoretical():
    colnames = ["Taxon","Theoretical Distribution","ds1_theoretical_reads","ds2_theoretical_reads"]
    tax1s = ["tax1 s", 0.2, 20000, 20000]
    tax2s = ["tax2 s", 0.1, 10000, 10000]
    theoretical_species = pd.DataFrame(columns = colnames, data=[tax1s,tax2s])
    tax1g = ["tax1", 0.5, 50000, 30000]
    tax2g = ["tax2", 0.2, 20000, 10000]
    theoretical_genus = pd.DataFrame(columns = colnames, data=[tax1g,tax2g])
    theoretical = {'species': theoretical_species, "genus": theoretical_genus}
    return theoretical


ident_krlines = [[0.2, 200, 200, 'S', '', 'tax1 s'],
                 [0.001, 1, 1, 'S', '', 'tax1 s p'], # diff strain of tax1 s
                 [0.1, 200, 200, 'S', '', 'tax2 s'],
                 [0.5, 501, 501, 'G', '', 'tax1'],
                 [0.1, 200, 200, 'G', '', 'tax2']]

zero_krlines = [[0, 0, 0, 'S', '', 'tax1 s'],
                 [0, 0, 0, 'S', '', 'tax1 s p'], # diff strain of tax1 s
                 [0, 0, 0, 'S', '', 'tax2 s'],
                 [0, 0, 0, 'G', '', 'tax1'],
                 [0, 0, 0, 'G', '', 'tax2']]

fp_krlines = [[0.2, 200, 200, 'S', '', 'tax3 s'],
                 [0.001, 1, 1, 'S', '', 'tax3 s p'], # diff strain of tax1 s
                 [0.1, 200, 200, 'S', '', 'tax4 s'],
                 [0.5, 201, 201, 'G', '', 'tax3'],
                 [0.1, 200, 200, 'G', '', 'tax4']]

def make_kreportD(krlines = ident_krlines, mult_by = None):
    if mult_by:
        new_kl = []
        for kl in krlines:
            mult_kl = [i * mult_by for i in kl[0:3]] + kl[3:]
            new_kl.append(mult_kl)
        krlines = new_kl
    colnames = ['proportion', 'cumulative_count', 'level_count', 'rank', 'taxid', 'name']
    kreportD = pd.DataFrame(columns = colnames, data=krlines)
    return kreportD


### TODO: I've just been visually checking these against my expectations.
# need to add these expectations into the tests so automatically checked!

def test_truepos_identical():
    theor = make_theoretical()
    tool1_kreport = make_kreportD()
    kr_dict = {"tool1": tool1_kreport, "Sourmash-k31": tool1_kreport}
    tr = []
    for rd in ["ds1"]:
        for tool, kr in kr_dict.items():
            taxR = assess_kreport("mc1", rd, theor, tool, kr, taxa_searchterms)
            tr+=taxR
    testD = pd.DataFrame.from_records(tr, columns=TaxResult._fields)
    print(testD)

def test_falseneg_zero():
    theor = make_theoretical()
    tool1_kreport = make_kreportD(krlines=zero_krlines)
    kr_dict = {"tool1": tool1_kreport, "Sourmash-k31": tool1_kreport}
    tr = []
    for rd in ["ds1"]:
        for tool, kr in kr_dict.items():
            taxR = assess_kreport("mc1", rd, theor, tool, kr, taxa_searchterms)
            tr+=taxR
    testD = pd.DataFrame.from_records(tr, columns=TaxResult._fields)
    print(testD)

def test_falseneg_detection_threshold():
    theor = make_theoretical()
    # multiply by 0.001 to drop below detection threshold (still TP for 'dany' column)
    tool1_kreport = make_kreportD(krlines=ident_krlines, mult_by = 0.001)
    kr_dict = {"tool1": tool1_kreport, "Sourmash-k31": tool1_kreport}
    tr = []
    for rd in ["ds1"]:
        for tool, kr in kr_dict.items():
            taxR = assess_kreport("mc1", rd, theor, tool, kr, taxa_searchterms)
            tr+=taxR
    testD = pd.DataFrame.from_records(tr, columns=TaxResult._fields)
    print(testD)

def test_all_falseneg_falsepos():
    theor = make_theoretical()
    tool1_kreport = make_kreportD(krlines=fp_krlines)
    #tool2_kreport = make_kreportD(mult_by=0.5)
    kr_dict = {"tool1": tool1_kreport, "Sourmash-k31": tool1_kreport}
    tr = []
    for rd in ["ds1"]:
        for tool, kr in kr_dict.items():
            taxR = assess_kreport("mc1", rd, theor, tool, kr, taxa_searchterms)
            tr+=taxR
    testD = pd.DataFrame.from_records(tr, columns=TaxResult._fields)
    print(testD)
