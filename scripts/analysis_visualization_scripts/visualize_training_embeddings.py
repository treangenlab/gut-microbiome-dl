import csv
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage
from matplotlib.patches import Patch

# create sets to contain the embeddings of each bioproject
bioprojects = dict({
    'PRJDB13214':set({}), 'PRJEB23010':set({}), 'PRJEB57404':set({}), 'PRJEB13870':set({}), 'PRJEB15257':set({}), 
    'PRJEB25008':set({}), 'PRJEB25727':set({}), 'PRJEB46960':set({}), 'PRJEB55125':set({}), 'PRJEB58436':set({}), 
    'PRJEB23147':set({}), 'PRJEB25514':set({}), 'PRJEB28701':set({}), 'PRJEB32135':set({}), 'PRJEB55713':set({}), 
    'PRJEB60573':set({}), 'PRJEB60773':set({}), 'PRJNA1049470':set({}), 'PRJDB11444':set({})
    })

# health categories
infectious_disease = set({'PRJDB13214', 'PRJEB23010', 'PRJEB57404'})
noninfectious_disease = set({'PRJEB13870', 'PRJEB15257', 'PRJEB25008', 'PRJEB25727', 'PRJEB46960', 'PRJEB55125', 'PRJEB58436'}) # PRJEB55125, PRJEB15257
healthy = set({'PRJEB23147', 'PRJEB25514', 'PRJEB28701', 'PRJEB32135', 'PRJEB55713', 'PRJEB60573', 'PRJEB60773', 'PRJNA1049470'})
other = set({'PRJDB11444'})


# sequencing method broad
novaseq = set({'PRJDB13214', 'PRJEB57404', 'PRJEB55125', 'PRJEB60573', 'PRJEB60773', 'PRJNA1049470', 'PRJDB11444'})
hiseq = set({'PRJEB23010', 'PRJEB13870', 'PRJEB25008', 'PRJEB25727', 'PRJEB46960', 'PRJEB58436', 'PRJEB23147', 'PRJEB25514', 'PRJEB28701', 'PRJEB55713'})
miseq = set({'PRJEB15257'})
nextseq = set({'PRJEB32135'})

# sequencing method specific
novaseq_6000 = set({'PRJDB13214', 'PRJEB57404', 'PRJEB55125', 'PRJEB60573', 'PRJEB60773', 'PRJNA1049470', 'PRJDB11444'})
hiseq_2000 = set({'PRJEB28701'})
hiseq_2500 = set({'PRJEB23010', 'PRJEB13870', 'PRJEB25727', 'PRJEB23147'})
hiseq_3000 = set({'PRJEB58436'})
hiseq_4000 = set({'PRJEB25008', 'PRJEB46960', 'PRJEB55713'})
hiseq_x_ten = set({'PRJEB25514'})
miseq = set({'PRJEB15257'})
nextseq_500 = set({'PRJEB32135'})


def read_embeddings_csv(fpath):
    """
    read the embeddings into a list
    """
    with open(fpath, newline='') as csvfile:
        reader = csv.reader(csvfile)
        next(reader)
        embeddings = []
        for row in reader:
            embeddings.append(tuple([float(item) for item in row]))
    return embeddings

def read_row_names(fpath):
    """
    read the file names that correspond to each embedding
    (in the same order) into a list
    """
    with open(fpath, 'r') as file:
        lines = file.readlines()
        return [line.strip().split('_')[0] for line in lines]

def sort_by_bioprojects(embeddings, row_bioprojects):
    """
    assign each embedding to its bioproject
    """
    for i in range(len(row_bioprojects)):
        bioprojects.get(row_bioprojects[i]).add(embeddings[i])

embeddings = read_embeddings_csv('path_to_embeddings.csv') #TODO: replace path
row_bioprojects = read_row_names('path_to_row_names.txt') #TODO: replace path -- row names in this file should not end with .txt
sort_by_bioprojects(embeddings, row_bioprojects)


def get_bioproject_category(bioproject):
    if bioproject in infectious_disease:
        return 'Infectious'
    elif bioproject in noninfectious_disease:
        return 'Noninfectious'
    elif bioproject in healthy:
        return 'Healthy'
    elif bioproject in other:
        return 'Other'
    else:
        return 'Unknown'
    
def get_sequencing_category(bioproject):
    if bioproject in novaseq:
        return 'NovaSeq'
    elif bioproject in hiseq:
        return 'HiSeq'
    elif bioproject in miseq:
        return 'MiSeq'
    elif bioproject in nextseq:
        return 'NextSeq'
    else:
        return 'Unknown'
    
def get_specific_sequencing_category(bioproject):
    if bioproject in novaseq_6000:
        return 'NovaSeq 6000'
    elif bioproject in hiseq_2000:
        return 'HiSeq 2000'
    elif bioproject in hiseq_2500:
        return 'HiSeq 2500'
    elif bioproject in hiseq_3000:
        return 'HiSeq 3000'
    elif bioproject in hiseq_4000:
        return 'HiSeq 4000'
    elif bioproject in hiseq_x_ten:
        return 'HiSeq X Ten'
    elif bioproject in miseq:
        return 'MiSeq'
    elif bioproject in nextseq:
        return 'NextSeq'
    else:
        return 'Unknown'


embeddings_in_label_order = []
labels = []
categories = []
sequencing_categories = []
specific_sequencing_categories = []

# for each embedding, categorize its bioproject. 
# these labels will be in the same order as the list of embeddings
for bioproject_id, embs in bioprojects.items():
    for emb in embs:
        embeddings_in_label_order.append(emb)
        labels.append(bioproject_id)
        categories.append(get_bioproject_category(bioproject_id))
        sequencing_categories.append(get_sequencing_category(bioproject_id))
        specific_sequencing_categories.append(get_specific_sequencing_category(bioproject_id))

# find cosine distance and cluster
X = np.array(embeddings_in_label_order)
cosine_dist = pdist(X, metric='cosine')
linkage_matrix = linkage(cosine_dist, method='average')
cosine_dist_matrix = squareform(cosine_dist)

# # -- BIOPROJECT --
# # assigns a color to each bioproject
# bioprojects_no_repeats = sorted(set(labels))
# palette = sns.color_palette("hls", len(bioprojects_no_repeats))
# project_to_color = dict(zip(bioprojects_no_repeats, palette))
# colors = [project_to_color[label] for label in labels]

# # plot it
# g = sns.clustermap(
#     cosine_dist_matrix,
#     row_linkage=linkage_matrix,
#     col_linkage=linkage_matrix,
#     cmap='viridis',
#     figsize=(10, 10),
#     col_colors=colors,
#     xticklabels=False,
#     yticklabels=False
# )

# legend_handles = [Patch(color=project_to_color[proj], label=proj) for proj in bioprojects_no_repeats]
# plt.suptitle("Samples Clustered and Colored by BioProject")
# plt.legend(handles=legend_handles, title='BioProject', bbox_to_anchor=(1.05, 1), loc='upper left')
# plt.show()


# -- HEALTH CATEGORIES --
# assigns a color to each bioproject
health_categories_no_repeats = sorted(set(categories))
palette = sns.color_palette("Set2", len(health_categories_no_repeats)) #set2 has nicer colors
category_to_color = dict(zip(health_categories_no_repeats, palette))
colors = [category_to_color[cat] for cat in categories]

# plot it
g = sns.clustermap(
    cosine_dist_matrix,
    row_linkage=linkage_matrix,
    col_linkage=linkage_matrix,
    cmap='viridis',
    figsize=(10, 10),
    col_colors=colors,
    xticklabels=False,
    yticklabels=False
)

legend_handles = [Patch(color=category_to_color[cat], label=cat) for cat in health_categories_no_repeats]
plt.suptitle("Samples Clustered and Colored by Health Category")
plt.legend(handles=legend_handles, title='Category', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.show()


# -- BROAD SEQUENCING METHOD --
broad_no_repeats = sorted(set(sequencing_categories))
palette = sns.color_palette("Set2", len(broad_no_repeats))
category_to_color = dict(zip(broad_no_repeats, palette))
colors = [category_to_color[cat] for cat in sequencing_categories]

# plot it
g = sns.clustermap(
    cosine_dist_matrix,
    row_linkage=linkage_matrix,
    col_linkage=linkage_matrix,
    cmap='viridis',
    figsize=(10, 10),
    col_colors=colors,
    xticklabels=False,
    yticklabels=False
)

legend_handles = [Patch(color=category_to_color[cat], label=cat) for cat in broad_no_repeats]
plt.suptitle("Samples Clustered and Colored by Sequencing Method")
plt.legend(handles=legend_handles, title='Category', bbox_to_anchor=(1, 1), loc='upper left')
plt.show()


# -- SPECIFIC SEQUENCING METHOD --
specific_no_repeats = sorted(set(specific_sequencing_categories))
palette = sns.color_palette("Set2", len(specific_no_repeats))
category_to_color = dict(zip(specific_no_repeats, palette))
colors = [category_to_color[cat] for cat in specific_sequencing_categories]

# plot it
g = sns.clustermap(
    cosine_dist_matrix,
    row_linkage=linkage_matrix,
    col_linkage=linkage_matrix,
    cmap='viridis',
    figsize=(10, 10),
    col_colors=colors,
    xticklabels=False,
    yticklabels=False
)

legend_handles = [Patch(color=category_to_color[cat], label=cat) for cat in specific_no_repeats]
plt.suptitle("Samples Clustered and Colored by Specific Sequencing Method")
plt.legend(handles=legend_handles, title='Category', bbox_to_anchor=(1, 1), loc='upper left')
plt.show()
