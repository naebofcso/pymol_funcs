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
    main_cluster =  '{}_{}'.format(obj_name, cluster_suffix)
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
    
    myspace = {'bfactor': [], 'Position':[]}
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
    b_df['Mean B'] = b_df[backup_struct_list].mean(axis = 1)
    b_df['Std. Dev. B'] = b_df[backup_struct_list].std(axis = 1)

    # Calculate normalized values to put it all on the same scale 
    # and compare differences more clearly
    for struct in backup_struct_list:
        norm_factor = protein_b_stats(stuct)[-1]
        b_df[struct + '_norm'] = b_df[struct]/norm_factor
    b_df['Norm. Mean B'] = b_df[
            [x + '_norm' for x in backup_struct_list]].mean(axis = 1)
    b_df['Norm. Std. Dev. B'] = b_df[
            [x + '_norm' for x in backup_struct_list]].std(axis = 1)
    
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

