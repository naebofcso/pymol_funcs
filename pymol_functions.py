'''
Author: Jacob Wolff


Description:

Useful functions related to pymol that can aid in calculations. Functionality
is not guaranteed because I am not a good programmer. In code comments should
help you understand what I'm trying to do so you can make things better as needed.
I used these in a Jupyter notebook running PyMol.

'''
# Drawing a mesh around Caver-generated spheres
# Requires that you identify the main cluster and prune it as needed
# Inputs:
# obj_name - the name of the object you ran caver on
# cluster_suffix - something like 't001_1'. This is the suffix made by Caver
# b - the sphere radius cutoff you are using for accetable spheres. Water is 1.4A
# cutoff - the distance in angstroms that atoms can be from the main cluster when
# you are adding to your final object


def draw_caver_mesh(obj_name, cluster_suffix, b=1.4, cutoff=3):

    main_cluster = '{}_{}'.format(obj_name, cluster_suffix)
    # select all atoms that have a B-factor greater than 1.4
    # (radius of water in angstroms)
    cmd.select('main_cluster', '{} and b > {}'.format(main_cluster, b))
    # for each remaining cluster create a new selection for it of atoms that
    # have b >1.4 and are within 3 angstroms of the main cluster
    cluster_list = ['main_cluster']
    # get the names of all clusters produced
    obj_list = cmd.get_object_list('({}_t0*_1)'.format(obj_name))
    # remove the main cluster so it doesn't get in the way
    obj_list.remove(main_cluster)
    # for each cluster, create new selection of spheres with radius > b
    # and within cutoff of the main cluster
    for cluster in obj_list:
        cmd.select('cluster_' + cluster[-1],
                   '({} w. {} of {}) and b > {}'.format(cluster, cutoff,
                                                        main_cluster, b))
        cluster_list.append('cluster_' + cluster[-1])
    # Create final selection
    cmd.select('final_selection', ' or '.join(cluster_list))
    cmd.hide('spheres')
    # Create new object
    cmd.create('final_object', 'final_selection')
    cmd.color('blue', 'final_object')
    # default is 2, higher quality smooths mesh more
    cmd.set('mesh_quality', '5')
    cmd.show('mesh', 'final_object')

# Function to get the average side chain B-factor. Useful if you want to look
# at an average B-factor for side chains. Requires that you have a structure
# loaded in PyMol.
# Input: a string representation of your PyMol selection
# Output: two lists. First list contains the position. Second list contains
# average B-factor for the side chain at that position


def get_avg_sc_b(selection):

    myspace = {'bfactor': [], 'Position': []}
    cmd.iterate(selection, 'bfactor.append(b)', space=myspace)
    cmd.iterate(selection, 'Position.append(resi)', space=myspace)

    # setup Loop
    avg_b = []
    sc_b = []
    sc_pos = []
    current_pos = myspace['Position'][0]

    # iterate through collected data
    for pos, b in zip(myspace['Position'], myspace['bfactor']):
        if pos != current_pos:
            avg_b.append(sum(sc_b)/len(sc_b))
            # append current position to create a list that doesn't have duplicates.
            sc_pos.append(pos)
            sc_b = []
            sc_b.append(b)
            current_pos = pos
        else:
            # If there is no change, simply "update" the position and append
            # that atom's b-factor
            current_pos = pos
            sc_b.append(b)

    return sc_pos, avg_b

# This function compares side chain b-factors and returns everything in a
# Pandas dataframe. Assumes you have imported pandas as 'pd'.
# Input:
# struct_list - list of strings containing object names for what you want to
# examine
# extra_select - string containing the selection modifiers you need to select
# the residues that you care about. I think it requires that you only select
# side chains (sc.). An example: ' and chain A and sc.'. The first space is
# important too.
# Currently the normalization is not good. It normalizes across your selection,
# not across your structure. Probably throws things off.
# Returns:
# - Pandas Dataframe containing the average side chain b-factor for the selection
# in each structure you input and columns containing Mean and Standard Deviation
# b-factors along with the normalized versions. Normalizes using the maximum
# b-factor from all amino acids
# Requires: get_avg_sc_b, protein_b_stats, Pandas


def bfac_comp(struct_list, extra_select):

    # Loop through all the structures to get the data and identify if there
    # is a "largest" structure
    backup_struct_list = struct_list.copy()
    bfac_dict = {}
    longest_count = 0
    longest_name = ''
    for struct in struct_list:
        # Relies on the get_avg_sc_b function to get the data
        position, bfac = get_avg_sc_b(struct + extra_select)
        bfac_dict[struct] = {x:y for x,y in zip(position, bfac)}
        if len(bfac) > longest_count:
            longest_name = struct
            longest_count = len(bfac)

    # Handle potential gaps. Untested!
    struct_list.remove(longest_name)
    for struct in struct_list:
        # use sets to speed up comparisons
        if set(bfac_dict[longest_name].keys()) == set(bfac_dict[struct].keys()):
            continue
        else:
            # Loop through the missing keys and fill in the value from the larger
            # structure
            missing_keys = set(bfac_dict[longest_name].keys()) - \
                set(bfac_dict[struct].keys())
            for key in missing_keys:
                bfac_dict[struct][key] = bfac_dict[longest_name][key]

    # set up the initial parts of the dataframe that aren't easily looped through
    b_df = pd.DataFrame({'Position': bfac_dict[longest_name].keys(),
                         longest_name: bfac_dict[longest_name].values()})
    # Fill the remaining data by looping
    for struct in struct_list:
        b_df[struct] = bfac_dict[struct].values()

    # Calculate some basic information
    b_df['Mean B'] = b_df[backup_struct_list].mean(axis=1)
    b_df['Std. Dev. B'] = b_df[backup_struct_list].std(axis=1)

    # Calculate normalized values to put it all on the same scale
    # and compare differences more clearly
    for struct in backup_struct_list:
        norm_factor = protein_b_stats(stuct)[-1]
        b_df[struct + '_norm'] = b_df[struct]/norm_factor
    b_df['Norm. Mean B'] = b_df[
            [x + '_norm' for x in backup_struct_list]].mean(axis=1)
    b_df['Norm. Std. Dev. B'] = b_df[
            [x + '_norm' for x in backup_struct_list]].std(axis=1)

    return b_df

# Function to get basic B-factor statistics for protein molecules using pandas 
# Input:
# - object name
# Output: pandas series produced by pandas.Series.describe()
# count
# mean
# std
# min
# 25%
# 50%
# 75%
# max


def protein_b_stats(obj_name):

    myspace = {'bfactor': []}
    # Get the b-factor for every protein atom
    cmd.iterate(obj_name + ' and polymer.protein', 'bfactor.append(b)',
                space=myspace)
    return pd.Series(myspace['bfactor']).describe()


'''
keyRes function is designed to get the original positions of amino acids in a
multiple sequence alignment given a a list of positions for the reference
sequence. This should help when you need to write about what the amino acid
position is in a homologous sequence.

**Warning!**
This function only works if your reference sequence is the first sequence in
the alignment.

Inputs:
    - A BioPython alignment object containing your MSA
    - A string with what your reference sequence id in it

output:
    - A dictionary where the keys are the sequence ids from the alignment object
    and the values are lists with the one letter code and the position such as
    'N100'. A dictionary was chosen to maximize output flexibility and avoid
    depending on external packages.
'''


def keyRes(alignment, refid: str, keyPos: list[int]) -> dict:

    # input sequence alignment object
    # Initiate a count for however many sequences are present
    num_seq = len(alignment[:, 0])
    counts = [0]*num_seq
    # get sequence ids from the alignment file: seq.id
    # use the sequence ids/names to make a dictionary
    aligned_pos = {x.id: [] for x in alignment}
    # Also keep a list of ids
    names = aligned_pos.keys()

    for pos in range(0, len(alignment[0])):
        # string representing the current column
        current_col = alignment[:, pos]
        # set False after going through all sequences
        importantPos = False
        for idx, res, seqname in zip(range(0, num_seq), current_col, names):
            nogap = res != '-'
            if nogap:
                counts[idx] += 1
            # assuming reference seq is first and gaps are not important
            if seqname == refid and counts[idx] in keyPos and nogap:
                importantPos = True
                aligned_pos[seqname].append(res+str(counts[idx]))
                continue
            if importantPos:
                if nogap:
                    aligned_pos[seqname].append(res+str(counts[idx]))
                else:
                    aligned_pos[seqname].append(res)
                continue

    return aligned_pos


'''
get_seq_offset is a function to handle differences in sequence between a structure
and its Uniprot sequence. It simply downloads the canonical PDB sequence and the
corresponding Uniprot sequence for the Uniprot entry and trys to use substring
matching to get the position using the string index() function. The core
assumption is that the reason the uniprot sequence is different from the PDB
sequence is that there is simply an offset. Thus finding the index of a
substring match will tell you what the offset is. If there isn't a match, the
substring will be shifted by 1 in case there is a deletion (insertions I'm not
sure about) and the index will be modified to reflect this so the offset
remains consistent. A further issue to watch out for is that PyMol can have
negative positions and 0 is included which makes any offset calculated off by
1.


This function was designed to only require the requests package installed to
minimize dependences.

Inputs:
    - A list of PDB ids for structures whose sequences you can to get offsets
    for.

output:
    - A dictionary where the keys are the PDB ids and the value is an int with
    the offset. A dictionary was chosen to maximize flexibility and this can
    easily be adapted to an ordered data structure as needed.
'''


def get_seq_offset(pdb_ids_list: list[str]) -> dict:

    import json
    import requests
    # Set base url
    rcsb_request_url = 'https://data.rcsb.org/graphql?query=<query>'
    # Open request template from file
    request_template = '{\
      polymer_entities(entity_ids:[<pdbid>]) {\
        rcsb_id\
        uniprots{\
            rcsb_uniprot_accession\
        }\
        entity_poly{\
            pdbx_seq_one_letter_code_can\
        }\
      }\
    }'
    # Format pbd id list for the query by capitalizing, adding "_1" for entity,
    # and encasing in double quotes
    pdbs = ['"' + x.upper() + '_1' + '"' for x in pdb_ids_list]
    # Join the list into a string
    request_payload = ','.join(pdbs)
    # Put this string into the request template
    new_request = request_template.replace('<pdbid>', request_payload)
    # Put request in URL and make request
    rcsb_request = requests.get(rcsb_request_url.replace('<query>', new_request))
    # Load Json response
    rcsb_content = json.loads(rcsb_request.text)
    # Extract Uniprot entry and canonical sequence from each pdbid.
    # Everything remains in order
    pdb_seqs = [x['entity_poly']['pdbx_seq_one_letter_code_can']
               for x in rcsb_content['data']['polymer_entities']]
    uniprot_ids = [x['uniprots'][0]['rcsb_uniprot_accession'][0]
                  for x in rcsb_content['data']['polymer_entities']]

    uniprot_url = 'https://rest.uniprot.org/uniprotkb/<entry>.fasta'
    seq_dict = dict()
    for pdbid, pdb_seq, uniprot_accession in zip(pdb_ids_list,
                                                 pdb_seqs, uniprot_ids):
        uniprot_request = requests.get(uniprot_url.replace('<entry>',
                                                           uniprot_accession))
        # Remove header a join together without any newline chars
        cut = uniprot_request.text.index('\n')
        current_seq = uniprot_request.text[cut:].replace('\n', '')

        # Scrub through to find using a window size of continuity
        continuity = 5
        match_found = False
        maximum_attempts = 10
        pdb_index = 0
        uniprot_index = 0
        if len(pdb_seq) <= len(current_seq):
            start_idx = 0
            attempts = 0
            while match_found is not True and attempts <= maximum_attempts:
                if pdb_seq[start_idx:start_idx + continuity] in current_seq:
                    pdb_index = current_seq.index(pdb_seq[
                        start_idx:start_idx + continuity])

                    match_found = True
                else:
                    start_idx += 1
                    attempts += 1

        elif len(pdb_seq) > len(current_seq):
            start_idx = 0
            attempts = 0
            while match_found is not True and attempts <= maximum_attempts:
                if current_seq[start_idx:start_idx + continuity] in pdb_seq:
                    uniprot_index = pdb_seq.index(current_seq[
                        start_idx:start_idx + continuity]) + start_idx

                    match_found = True
                else:
                    start_idx += 1
                    attempts += 1

        else:
            seq_dict[pdbid] = 'check resseq value'
        if pdb_index != 0:
            seq_dict[pdbid] = pdb_index
        elif uniprot_index != 0:
            seq_dict[pdbid] = uniprot_index*-1
        else:
            seq_dict[pdbid] = 0

    return seq_dict


'''
msa_struct_compare is a function to assess whether your sequence alignments
make sense structurally. This is done by supplying positions for each protein,
reference positions, a reference structure name, and a cutoff. Assuming that
your structures aligned already, a dictionary containing the distance
comparisons will be output and distances larger than the cutoff will be printed
for immediate assessment.

This is the follow-up function to keyRes and get_seq_offset.

Inputs:
    - Pandas DataFrame containing the amino acid labeling information. Assumes
      the index contains the protein names and that the columns are each
      similar position. Assumes that the amino acids are in the form 'M100' for
      example.
    - List of query chains if you need to compare different chains per query
      structure
    - List of reference positions also in the form 'M100'
    - Reference structure name
    - Reference chain name
    - Distance cutoff in angstroms

Output:
    - Dictionary containing the comparisons. The key is the hyphenated
      comparison 'Query-reference'. Values are lists of tuples containing the
      query position, reference position, and average distance between the
      atoms.
'''


def msa_struct_compare(key_sites_df: pd.DataFrame, query_chains: list[str] = [],
                       ref_positions: list[str] = [], ref_struct: str = '',
                       ref_chain: str = '', cutoff: int = 5) -> dict:

    # Create dictionary to hold values
    distance_dict = {}

    # Roughly handle inputs that are lacking
    if ref_chain == '':
        ref_chain = 'A'
    if len(ref_positions) < 1:
        raise ValueError('No reference positions provided!')
    if len(query_chains) < 1:
        query_chains = ['A']*len(key_sites_df.index)
    if ref_struct == '':
        raise ValueError('No reference structure provided!')


    pdb_renumbering = [cmd.get_model(x).atom[0].resi for x in key_sites_df.index]
    for imp_res, ref_res in zip(key_sites_df.columns, ref_positions):
        for pos, query_struct, query_chain in zip(
                key_sites_df[imp_res], key_sites_df.index, query_chains):
            comparison_name = query_struct + '-' + ref_struct
            if comparison_name not in distance_dict:
                distance_dict[comparison_name] = []
            # cmd.distance returns the average distance across atoms
            average_dist = cmd.distance(name=comparison_name,
                            selection1='{} and chain {} and resi {}'.format(
                                        ref_struct, ref_chain, ref_res[1:]),
                            selection2='{} and chain {} and resi {}'.format(
                                        query_struct, query_chain, pos[1:]))
            if average_dist > cutoff:
                print('{} to {} distance betwen {} and {} is {}'.format(
                    ref_struct, query_struct, pos, ref_res, average_dist))
            distance_dict[comparison_name].append((pos, ref_res, average_dist))
            # Delete distance object that's created
            cmd.delete(comparison_name)
    return distance_dict

'''
match_struct_seq is a function to match protein structure sequence numbering to
their respective Uniprot sequence numbering to make selecting residues
standardized to the uniprot numbering. Requires that you have your structures
already loaded in pymol. Also requires data from get_seq_offset and keyRes.

Input:
    - List of protein chains for getting the sequence
    - A dictionary for conversions. My previous functions take in PDB ids but I
      load the structures into PyMol and rename to the protein names. This
      required that I convert back from protein name to PDB id to access the
      data from those functions.
    - A dictionary with offsets for the proteins that you are examining
    - Boolean indicating whether you need to convert or not. Default is True
      since this work is based on my previous functions

Output:
    None. This function simply modifies the residue numbering in the loaded
    structures based on the offset.
'''

def match_struct_seq(chains_list: list[str], conversion_dict: dict,
                     offset_dict: dict, convert: bool = True):

    for obj, chain in zip(cmd.get_object_list(), chains_list):
        offset = 0
        if convert is True:
            converted_id = conversion_dict[obj]
            offset = offset_dict[converted_id]
        else:
            offset = offset_dict[obj]

        myspace = {'resilist': []}
        cmd.iterate('{} and chain {} and n. CA'.format(obj, chain),
                    'resilist.append(resi)',
                    space=myspace)
        seq_idx = int(myspace['resilist'][0])
        if seq_idx == 1 and offset != 0:
            cmd.alter(obj + ' and polymer', 'resv += {}'.format(offset))
        elif abs(seq_idx - offset) > 0:
            if offset != 0:
                cmd.alter(obj + ' and polymer', 'resv += {}'.format(
                    abs(seq_idx - offset)-1))
            elif seq_idx < 0:
                cmd.alter(obj + ' and polymer', 'resv += {}'.format(
                    abs(seq_idx - offset)+1))
            elif seq_idx > 0:
                continue
