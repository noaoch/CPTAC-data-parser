import numpy
import pyranges as pr
import pybedtools
import pandas as pd
import cptac
from typing import Union
import copy
import logging
import mygene


############################################################################################

def parse_start_end(s: str) -> (int, int):
    """
    Given a string in the format ">chrX:start-end" from pybedtools, where start and end are in 0-coordinate system,
    extract the start and end positions, and record them as a tuple (start, end).

    :param s: string with format >chrX:start-end in 0-coordinate system
    :return: int tuple of 2 elements in format (start, end) in 0-coordinate system
    """
    s = s.split(":")
    assert(len(s) == 2)

    # Focus on the part after chrX
    s = s[1]

    # split to start and end positions
    s = s.split("-")
    assert(len(s) == 2)

    start = int(s[0])
    end = int(s[1])

    return start, end




##########################################################################################################################################################
##########################################################################################################################################################
##########################################################################################################################################################
##########################################################################################################################################################
##########################################################################################################################################################




######################## Sequence grouper for patients ##############################

def group_by_sample_id(maf: pd.DataFrame):
    """
    Given a maf dataframe, group the dataframe by sample id, and return a dictionary where the key is sample id,
    :param maf: mutation maf
    :return: dictionary where the key is sample id, and the value is the index of the maf dataframe

    >>> cptac.download("colon")
    >>> c = cptac.Colon()
    >>> maf = c._maf
    >>> result = group_by_sample_id(maf)
    >>> type(result["11CO059"])
    >>> len(result)

    >>> cptac.download("endometrial")
    >>> en = cptac.Endometrial()
    >>> en_maf = en._maf
    >>> result = group_by_sample_id(en_maf)
    >>> result
    """
    try:
        grouped = maf.groupby("SampleID")
    except:
        grouped = maf.groupby("Tumor_Sample_Barcode")
    result = {}
    for name, group in grouped:
        try:
            group = group.drop_duplicates(subset=['Start', 'Gene'])
            group = group.drop_duplicates(subset=['Protein_Change', 'Gene']) # drop the mutations that would result in the same protein change
            result[name] = group
        except: # files differ by their format slightly (in endometrial start is named Start_Position)
            try:
                group = group.drop_duplicates(subset=['Start_Position', 'Gene'])
                result[name] = group
            except:
                group = group.drop_duplicates(subset=['Start_Position', 'Hugo_Symbol'])
                result[name] = group

    return result


###############################   Revised One_mutation introduction   ########################

def one_mutation(ref1, ref2, patient, seen) -> Union[str, int, None]:
    """
    introduce one mutation
    >>> cptac.download("colon")
    >>> c = cptac.Colon()
    >>> maf = c._maf
    >>> p_zero = maf.iloc[0]
    >>> hg38 = read_hg38_v36()
    >>> hg19 = read_hg19()
    >>> s = {}
    >>> mutated_seq = one_mutation(hg38, hg19, p_zero, s)
    >>> original_gene = fetch_gene(hg19, 'A1BG', 'chr19')
    >>> original_seq = gene2seq(original_gene[0], 'hg19')
    >>> start, end = parse_start_end(original_seq[0])
    >>> original_seq = original_seq[1]
    >>> mut_index = (58864357 - 1) - start
    >>> mutated_seq[mut_index - 3:mut_index + 3]
    >>> original_seq[mut_index - 3:mut_index + 3]

    >>> maf_dic = group_by_sample_id(maf)
    >>> p_one = maf_dic["01CO014"]
    >>> p_one = p_one[(p_one['Gene'] == 'CFP')].iloc[1]
    >>> print(f"before:{p_one['Ref']}")
    >>> print(f"after:{p_one['Alt']}")
    >>> mutated_seq_2 = one_mutation(hg38, hg19, p_one, {})
    >>> mutated_seq_2[1905:1915]

    >>> p_two = maf_dic["05CO037"]
    >>> p_two = p_two[(p_two['Gene'] == 'NDST2') & (p_two['Ref'] == 'AAAAAAA') & (p_two['Alt'] == 'AAAAAA')].iloc[0]
    >>> p_two
    >>> print(f"before:{p_two['Ref']}")
    >>> print(f"after:{p_two['Alt']}")
    >>> mutated_seq_2 = one_mutation(hg38, hg19, p_two, {})
    >>> mutated_seq_2
    """

    info = trial_reader(patient)
    symbol, build, mut_start, mut_end, before, after, chr = info
    ref = ref1 if build in ('GRCh38', 'hg38') else ref2 # Use hg38 or hg19

    locations_dict = {}  # New
    if symbol in seen:
        original, original_start, original_end, seq_to_mut, offset = seen[symbol]
        locations_dict[0] = (original, original_start, original_end, seq_to_mut, offset)  # New


    else:
        offset = {} # offset is a dict {mut end location in original seq: how many letters deviated from original seq after here}
        genes_to_mut = fetch_gene(ref, symbol, chr)  # This is a list of gene locations to try out
        if genes_to_mut is None:
            return None

        # Get the sequence for every possible gene location
        for i in range(len(genes_to_mut)):  # New
            seq_to_mut = gene2seq(genes_to_mut[i], build)  # Fetch the sequence for the ith gene location
            original_start, original_end = parse_start_end(seq_to_mut[0])
            seq_to_mut = seq_to_mut[1]  # This is the original seq since this is the first time we saw this gene id
            original = seq_to_mut
            locations_dict[i] = (original, original_start, original_end, seq_to_mut, offset)


    for key in locations_dict:
        original_seq, original_start, original_end, seq_to_mut, offset = locations_dict[key]
        print(f"offset before muattion: {offset}")
        indices = get_start_end_index(original_start, original_end, mut_start, mut_end)
        if indices is None:
            continue  # Try the next possible location
        start_index, end_index = indices

        x = mutate(seq_to_mut, offset, start_index, end_index, before, after)
        if x is None:
            continue  # Try next location
        mutated_seq, new_offset = x
        print(f"offset after mutation: {new_offset}")
        print(f"start index is {start_index}")
        print(f"end index is {end_index}")

        # Update seen
        seen[symbol] = original_seq, original_start, original_end, mutated_seq, new_offset
        return mutated_seq

    return None  # Return None if no possible location matches



    # print(f"offset before muattion: {offset}")
    # indices = get_start_end_index(original_start, original_end, mut_start, mut_end)
    # if indices is None:
    #     return None
    # start_index, end_index = indices
    #
    # x = mutate(seq_to_mut, offset, start_index, end_index, before, after)
    # if x is None:
    #     return -1
    # mutated_seq, new_offset = x
    # print(f"offset after mutation: {new_offset}")
    # print(f"start index is {start_index}")
    # print(f"end index is {end_index}")
    #
    # # Update seen
    # seen[symbol] = original, original_start, original_end, mutated_seq, new_offset
    # return mutated_seq



def get_start_end_index(gene_start, gene_end, mut_start, mut_end):
    """
    Get the start and end indices in the gene sequence where mutation occurred (as a tuple), where gene_start and
    gene_end are based on the original, un-mutated fresh gene.
    """
    try:
        assert(gene_start <= mut_start <= mut_end <= gene_end)
    except AssertionError:
        print(f"gene_start: {gene_start}, mut_start: {mut_start}, mut_end: {mut_end}, gene_end: {gene_end}")
        print("This sequence has mutation outside of the reference segment.")
        return None

    start_index = mut_start - gene_start
    end_index = mut_end - gene_start
    return start_index, end_index


############# Trial Reader #############
def trial_reader(patient: pd.Series) -> tuple[str, str, int, int, str, str, str]:
    """
    Read the patient file and output gene_id, build, mutation_start, mutation_end, before, after, chr as a tuple
    NOTE: start positions are converted to 0-based coordinates HERE!!!
    :param patient:
    :return:

    >>> cptac.download("colon")
    >>> c = cptac.Colon()
    >>> maf = c._maf
    >>> p_zero = maf.iloc[0]
    >>> gr = pr.read_gtf("gencode.v22.annotation.gtf")
    >>> v22 = gr.df
    >>> v22["gene_id"] = v22["gene_id"].apply(lambda x: x.split(".")[0])

    """
    # GENE_SYMBOL
    # This gene symbol is geared towards the colon dataset
    try:
        gene_symbol = patient["Gene"]
    except:
        # Hugo_symbol is first seen in lscc
        gene_symbol = patient["Hugo_Symbol"]

    if gene_symbol.startswith("ENSG"):  # if this is not the gene symbol, but the gene id
        # This column notation is first found in endometrial
        gene_symbol = patient["SYMBOL"]
    print(f"Gene Symbol: {gene_symbol}")

    try: # BUILD
        ver = patient["NCBI_Build"]
    except KeyError: # If build not found
        ver = "GRCh37" # Assume we are using hg19 (First case where the build is not provided is in colon)
    print(f"Build: {ver}")


    # Mutations all use 1-based coordinates, convert them to 0-based
    try:  # START
        mut_start = int(patient["Start_Position"]) - 1
    except KeyError:
        # this notation first observed in colon
        mut_start = int(patient["Start"]) - 1

    try:  # END
        mut_end = int(patient["End_Position"])
    except KeyError:
        try:
            # This notation is first observed in lscc
            mut_end = int(patient["End_position"])
        except KeyError:
            # This notation is first observed in colon
            mut_end = int(patient["End"])

    try: #BEFORE-AFTER
        mut_before = patient["Reference_Allele"]
        mut_after = patient["Tumor_Seq_Allele2"]
    except KeyError:
        # First observed in colon
        mut_before = patient["Ref"]
        mut_after = patient["Alt"]


    try: # Find on which chromosome the mutation occurred
        chr = patient["Chromosome"]
    except KeyError:
        # In the colon file chr is recorded as indices in the dataframe, access it
        chr = str(patient.name)

    # Some files only state the chromosome number X, while we want chrX
    if not chr.startswith("chr"):
        chr = "chr" + chr

    return gene_symbol, ver, mut_start, mut_end, mut_before, mut_after, chr







############# Revised Sequence Fetcher ##############
def fetch_gene(ref, gene_symbol, chr):
    """
    Return a dataframe for a gene in ref w/ specified gene_id and chr
    """
    genes = []

    whole_seq = ref[(ref["gene_name"] == gene_symbol) & (ref["Chromosome"] == chr) & (ref["Feature"] == "gene")]
    print(f"chromosome given: {chr}")
    if whole_seq.empty: # See if whole_seq is empty (i.e. the gene_id is not found in ref)
        print(f"gene symbol {gene_symbol} not found in ref")
        return None

    # we could see multiple genes with the feature 'gene',
    # try all of them
    for j in range(len(whole_seq)):
        jth_row = whole_seq.iloc[j].to_frame().transpose()
        cleaned_gene = jth_row.loc[:, ["Chromosome", "Start", "End"]]

        if type(cleaned_gene['Start'].iloc[0]) == numpy.int64:  # if hg38
            new_gene = copy.deepcopy(cleaned_gene)
            genes.append(new_gene)

        else:
            # Start is a tuple for hg19, so we have multiple genes
            print(f"start: {cleaned_gene['Start'].iloc[0]}")
            for i in range(len(cleaned_gene['Start'].iloc[0])):
                new_gene = copy.deepcopy(cleaned_gene)
                try:
                    new_gene['Start'] = new_gene['Start'].iloc[0][i]
                    new_gene['End'] = new_gene['End'].iloc[0][i]
                except:
                    print("Error in fetch_gene")
                    print(new_gene)
                    print(new_gene['Start'])
                    print(new_gene['End'])
                    print(i)
                    print(len(cleaned_gene['Start']))
                    raise
                # The Start and End here are 1-based, but pybedtools uses 0-based coordinates here, so change them
                new_gene['Start'] = new_gene['Start'].apply(lambda x: x - 1)
                # No need to change End, because we want to include the last base pair of the gene segment
                genes.append(new_gene)

    return genes






################## Revised gene -> seq ###############
def gene2seq(gene, build) -> tuple[str, str]:
    """
    Convert the given gene with entries: "chromosome, start, end", to the actual sequence
    """
    s = gene.to_string(header=False, index=False, index_names=False)
    a = pybedtools.BedTool(s, from_string=True)

    # Assign sequence
    a = a.sequence(fi="GRCh38.d1.vd1.fa") if build in ['GRCh38', 'hg38'] else a.sequence(fi="GENCODE fasta and GTF/GRCh37.p13.genome.fa")
    desired_seq = (open(a.seqfn).read())
    split_seq = desired_seq.split('\n')

    # There's a newline at the end of desired_seq, so the last element of split_seq is an empty string (remove it!)
    split_seq = split_seq[:-1]
    assert (len(split_seq) == 2)
    return split_seq[0], split_seq[1] # Format: (chrX:Start-End, actual seq)




################# Revised mutate letter-seq ############
def mutate(seq_to_mut, offset: dict[int, int], start_index, end_index, before, after):
    """
    We get a seq to mutate, the current offset from original gene, start&end indices, before&after.
    Return a tuple (mutated sequence, new offset)

    >>> s1 = 'TTTTTTTT'
    >>> z, w = mutate(s1, {}, 0, 1, 'T', '-')
    >>> z
    >>> w
    >>> z, w = mutate(z, w, 7, 8, 'T', '-')
    >>> z
    >>> w
    >>> z, w = mutate(z, w, 2, 3, 'T', '-')
    >>> z
    >>> w  # ERROR: the higher indices are compounded from before, they forgot what they originally were

    >>> s2 = 'AAAATTTT'
    >>> z, w = mutate(s2, {}, 0, 1, '-', 'C')
    >>> z
    'CAAAATTTT'
    >>> w
    {1: 1}
    >>> z, w = mutate(z, w, 1, 7, 'AAATTT', 'AAATTTC')
    >>> z
    'CAAAATTTCT'
    >>> w
    {1: 1, 7: 2}

    """
    # offset is sorted dict, go until we see the most appropriate (tightest lowerbound) end location
    cur_best = 0 # start with no offset from original seq (i.e. so far identical to original seq)
    for k in offset:
        if k <= start_index: # if a mutation ended right before this one
            cur_best = offset[k] # update the best numerical offset so far
        else:
            break
    # print(f"best:{cur_best}")

    adjusted_start_index = cur_best + start_index # seq_to_mut may not be the original gene, yet start&end indices are based on original gene, so make adjustment
    adjusted_end_index = cur_best + end_index

    if before == '-':
        # mutation is an insertion
        mutated_seq = seq_to_mut[:adjusted_start_index] + after + seq_to_mut[adjusted_start_index:] # for insertion, start == end in the original file
        mutated_offset = len(after) + cur_best # the newly mutated gene is longer than what we had before
        offset[end_index] = mutated_offset # end index in original seq leads to a new offset threshold
        delta = len(after)  # record change in length starting from this point (for updating mutation checkpoints in the future)

    elif after == '-':
        # mutation is a deletion
        try:
            assert (seq_to_mut[adjusted_start_index:adjusted_end_index] == before)
        except AssertionError:
            print("This sequence has mutation unmatched to the reference sequence")
            print(f"seq neighborhood: {seq_to_mut[adjusted_start_index-4:adjusted_end_index+5]}")
            print(f"mutation: {before} -> {after}")
            # exit(-1)
            return None
        mutated_seq = seq_to_mut[:adjusted_start_index] + seq_to_mut[adjusted_end_index:]
        mutated_offset = cur_best - len(before) # newly mutated seq shorter than before
        offset[end_index] = mutated_offset
        delta = 0 - len(before)

    else:
        # mutation is a substitution
        try:
            assert (seq_to_mut[adjusted_start_index:adjusted_end_index] == before)
        except AssertionError:
            print("This sequence has mutation unmatched to the reference sequence")
            print(f"problematic seq: {seq_to_mut[adjusted_start_index:adjusted_end_index]}")
            print(f"seq neighborhood: {seq_to_mut[adjusted_start_index-4:adjusted_end_index+4]}")
            print(f"mutation: {before} -> {after}")
            # exit(-1)
            return None
        mutated_seq = seq_to_mut[:adjusted_start_index] + after + seq_to_mut[adjusted_end_index:]
        mutated_offset = (len(after) - len(before)) + cur_best  # seq turned shorter -> neg offset, otherwise positive
        offset[end_index] = mutated_offset
        delta = len(after) - len(before)

    sorted_keys = list(offset.keys()) # attempt to sort the offset dict
    sorted_keys.sort()
    # Every offset after this offset is also effected
    sorted_offset = {}
    for k in sorted_keys:
        v = offset[k] + delta if k > end_index else offset[k]
        sorted_offset[k] = v

    return mutated_seq, sorted_offset







###############################  Mutate sequences for all patients (non-colon)  ###################################

def introduce_all_mutation(ref1: pd.DataFrame, ref2: pd.DataFrame, patients_maf: pd.DataFrame) -> pd.DataFrame:
    """
    For each patient mutation in patients_maf, introduce it to reference sequence, and keep only the exon/UTR regions.
    Return a dataframe where each column is a patient, and each box on that column is an exon/UTR region.

    :param ref1: reference sequence for GRCh38
    :param ref2: reference sequence for GRCh37
    :param patients_maf: dataframe recording all patients (for a particular mutation such as endometrial)
    :return: dataframe where each column is a patient, and each box on that column is an exon/UTR region.
    """
    # ref: gene id -> sequence
    # seen_gene: gene_symbol -> (original, original_start, original_end, currently mutated, offset)
    j = 0
    d = {}
    failed_genes = 0
    grouped = group_by_sample_id(patients_maf)

    for key in grouped: # We are working on the same patient
        seen_gene = {}
        patient_df = grouped[key]
        for i in range(len(patient_df)):
            ith_mutation = patient_df.iloc[i]
            mutated_seq = one_mutation(ref1, ref2, ith_mutation, seen_gene)

            ########################## Error checking
            if mutated_seq is None:
                failed_genes += 1
            ##########################
            else:
                print(f"patient key: {key}")
                print("Done!")

        patient_j_dic = {}  # Make 1 column for each patient
        for k in seen_gene:
            patient_j_dic[k] = seen_gene[k][3]
        d[f"patient{j}"] = patient_j_dic
        print(f"####################################### finished patient with ID: {key} ###########################################")
        j += 1

    print(f"failures: {failed_genes}")
    return pd.DataFrame(d)






def focus_on_gene_name(gtf: pd.DataFrame) -> pd.DataFrame:
    """
    Focus on the gene_name column by renaming it to "gene_id", and remove the gene_id column.
    :param gtf: the gtf dataframe
    :return: the gtf dataframe with gene_name column renamed to "gene_id", and gene_id column removed.
    """
    gtf = gtf.drop(columns=["gene_id"])
    gtf = gtf.rename(columns={"gene_name": "gene_id"})

    return gtf


def read_hg19() -> pd.DataFrame:
    """
    In the legacy hg19 file, the "Gene" column is written in the formal "HGNC_Gene_symbol|Gene_id"
    Parse the column to contain only the Gene Symbol

    >>> x = read_hg19()
    >>> x
    >>> type(x.loc[0, 'Start'])
    >>> x['Start'] = x['Start'].apply(lambda y: int(y))
    >>> x
    """
    # TODO: Parse "Gene" column to contain only HGNC Gene Symbol.
    #  Add columns Start, End, and Chromosome based on column "GeneLocus".
    #  Rename column "FeatureType" to "Feature"

    hg19 = pd.read_csv('TCGA.hg19.June2011.gaf', sep='\t')
    rows_to_keep = [x for x in range(len(hg19)) if type(hg19.iloc[x]['GeneLocus']) != float]  # Drop unwanted rows
    hg19 = hg19.iloc[rows_to_keep, :]
    hg19['Gene'] = hg19['Gene'].apply(lambda x: str(x).split('|')[0])
    hg19['Chromosome'] = hg19['GeneLocus'].apply(lambda x: str(x).split(':')[0])

    hg19['Start'] = hg19['GeneLocus'].apply(lambda x: map_start(x))
    hg19['End'] = hg19['GeneLocus'].apply(lambda x: map_end(x))
    # hg19 = hg19[(len(hg19['Start']) != 0)]
    #
    # hg19['Start'] = hg19['Start'].apply(lambda x: int(x))
    # hg19['End'] = hg19['End'].apply(lambda x: int(x))
    hg19 = hg19.rename(columns={"FeatureType": "Feature"})
    hg19 = hg19.rename(columns={"Gene": "gene_name"})
    hg19.drop_duplicates(subset=['Feature', 'Start', 'End'], inplace=True)

    return hg19


def map_start(s: str) -> tuple[int]:
    x = s.split(';')
    lst = []
    for element in x:
        if len(str(element).split(':')) > 1:
            start_loc = str(element).split(':')[1].split('-')[0]
            lst.append(int(start_loc))
    return tuple(lst)

def map_end(s: str) -> tuple[int]:
    x = s.split(';')
    lst = []
    for element in x:
        if len(str(element).split(':')) > 1:
            end_loc = str(element).split(':')[1].split('-')[1]
            lst.append(int(end_loc))
    return tuple(lst)






def read_hg38_v36() -> pd.DataFrame:
    """
    Return the annotation dataframe for hg38 v36 (obtained from GDC)
    :return:

    >>> x = read_hg38_v36()
    >>> x
    """
    gr = pr.read_gtf("gencode.v36.annotation.gtf")
    hg38 = gr.df
    # In the Ensemble gene id, we have the format "ENSG00000xxxxxx.y", where y means the version number it is on.
    # We want to remove the version number, so we split the string by ".", and take the first element.
    hg38["gene_id"] = hg38["gene_id"].apply(lambda x: x.split(".")[0]) # <- Is this line really necessary???
    return hg38




# def drop_garbage(maf_dic):
#     """
#     Drop contradictory mutations from dic (the ones where their mutation coordinates intersect)
#
#     >>> cptac.download("colon")
#     >>> colon = cptac.Colon()
#     >>> colon_maf = colon._maf
#
#     """
#     for key in maf_dic:
#         p = maf_dic[key]
#         bad = pd.DataFrame(columns=p.columns)
#         for i in range(len(p)):
#             ith_entry = p.iloc[i]
#             for j in range(i + 1, len(p)):
#                 jth_entry = p.iloc[j]
#                 if (ith_entry['Start'] <= jth_entry['Start'] and ith_entry['End'] >= jth_entry['End']) or \
#                         (ith_entry['Start'] >= jth_entry['Start'] and ith_entry['End'] <= jth_entry['End']):
#                     bad = pd.concat([bad, ith_entry])
#                     bad = pd.concat([bad, jth_entry])
#
#         p = p[~p.x1.isin(bad.x1)]
#         maf_dic[key] = p
#     return maf_dic







if __name__ == "__main__":
    # import doctest
    # doctest.testmod()

    # colon_maf = colon_maf[(colon_maf['Gene'] == 'SPAG11B')]
    # # colon_maf = colon_maf[(colon_maf['Ref'] == 'A') & (colon_maf['Alt'] == 'G')]
    # colon_maf = colon_maf[(colon_maf['SampleID'] == '05CO041')]
    #
    # colon_maf = colon_maf.drop_duplicates(subset=['Start', 'Gene'])
    # colon_maf = colon_maf.drop_duplicates(subset=['Protein_Change', 'Gene'])  # drop the mutations that would result in the same protein change
    # assert 0

    print("Initializing")
    hg38 = read_hg38_v36()
    hg19 = read_hg19()

    # # testing sequence builders for exons/UTRs, and whole sequence
    # print("building dict...")
    # x = annotations_df.loc[44, "gene_id"]
    # x = x.split(".")[0]
    # d = build_exonUTR_sequence_string(annotations_df, x)
    # print(d)
    # #
    # print("building whole sequence string...")
    # c = build_whole_sequence_string(annotations_df, x)
    # print(c)
    #
    # ref_start_end = parse_start_end(c[0])
    # ref_start = ref_start_end[0]
    # ref_end = ref_start_end[1]
    #
    # for key in d:
    #     exon_utr_start_end = parse_start_end(key)
    #     exon_utr_start = exon_utr_start_end[0]
    #     exon_utr_end = exon_utr_start_end[1]
    #
    #     start_index = exon_utr_start - ref_start
    #     end_index = exon_utr_end - ref_start
    #
    #     # Check the exon / utr can be correctly indexed on the reference sequence
    #     # Note that the end index is not inclusive
    #     # (i.e. our exon / UTR ends with the base immediately before the end index)
    #     print(c[1][start_index:end_index])
    #     print(d[key])
    #     assert (c[1][start_index:end_index] == d[key])

    # # testing parse_start_end
    # keys = list(d.keys())
    # k = keys[0]
    #
    # tup = parse_start_end(k)
    # print(tup)
    #
    # print(d[k])
    # print(len(d[k]))




    print("using cptac...")
    # cptac.download("endometrial")
    # en = cptac.Endometrial()
    # en_maf = en._maf

    # cptac.download("brca")
    # brca = cptac.Brca()
    # brca_maf = brca._maf

    # cptac.download("ccrcc")
    # ccrcc = cptac.Ccrcc()
    # ccrcc_maf = ccrcc._maf

    # cptac.download("gbm")
    # gbm = cptac.Gbm()
    # gbm_maf = gbm._maf

    # cptac.download("hnscc")
    # hnscc = cptac.Hnscc()
    # hnscc_maf = hnscc._maf

    # cptac.download("lscc")
    # lscc = cptac.Lscc()
    # lscc_maf = lscc._maf

    # cptac.download("luad")
    # luad = cptac.Luad()
    # luad_maf = luad._maf

    # cptac.download("ovarian")
    # ov = cptac.Ovarian()
    # ov_maf = ov._maf

    cptac.download("pdac")
    pdac = cptac.Pdac()
    pdac_maf = pdac._maf

    # cptac.download("colon")
    # colon = cptac.Colon()
    # colon_maf = colon._maf

    # # check if there's anything in the negative strand (Not suitable for colon)
    # print("checking negative strand...")
    # colon_maf_neg = colon_maf[colon_maf["Strand"] == "-"]
    # print(colon_maf_neg)
    #
    # assert colon_maf_neg.empty

    # Parse the data
    print("introducing all mutations...")
    pdac_mutations_df = introduce_all_mutation(hg38, hg19, pdac_maf)
    pdac_mutations_df.to_csv("pdac_mutations.csv")
    # df = pd.read_csv("brca_mutations.csv")
























    # # testing introduce_one_mutation

    # debug log:
    # 1. Patient 381: mutation gene id is not found in the reference sequence
    # 2. Patient 977: mutation is outside the reference segment
    # 3. Patient 1466: mutation has gene_id that matches 2 genes (different chromosomes) in the reference sequence

    # patient_zero = brca_mutations.iloc[58]
    # mutated_seq = introduce_one_mutation(hg38, hg19, patient_zero)

    # Take note that pybedtools uses 0-based coordinate system, while gtf and maf use 1-based coordinate system


