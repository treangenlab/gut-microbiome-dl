import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
import seaborn as sns
from scipy.stats import ttest_ind


df = pd.read_csv('path_to_diabimmune_metadata_and_embeddings_merged.csv') #TODO: replace path

# find the correlation between node 27 and node 44's activation values 
# and participant's age at collection. Using spearman's rank correlation 
# coefficient (over something like Pearson's) because we want to 
# look at monotonic (but not strictly linear) relationships

# node 27's correlation with age
corr, p_val = spearmanr(df['27'], df['Age_at_Collection'])
print("node 27 with age stats:")
print(corr)
print(p_val)

plt.figure(figsize=(8, 6))
plt.figure(figsize=(8, 6))
sns.scatterplot(
    x='Age_at_Collection',
    y='27',
    data=df,
    s=60,
    alpha=0.7
)
plt.title(f"Age and Node 27 Value Relationship")
plt.xlabel("Age at Collection", fontsize=12)
plt.ylabel("Value of Node 27", fontsize=12)
plt.grid(True)
plt.tight_layout() # looks strange otherwise
plt.show()


# node 44's correlation with age
corr, p_val = spearmanr(df['44'], df['Age_at_Collection'])
print("node 44 with age stats:")
print(corr)
print(p_val)

plt.figure(figsize=(8, 6))
sns.scatterplot(
    x='Age_at_Collection',
    y='44',
    data=df,
    s=60,
    alpha=0.7
)
plt.title(f"Age and Node 44 Value Relationship")
plt.xlabel("Age at Collection", fontsize=12)
plt.ylabel("Value of Node 44", fontsize=12)
plt.grid(True)
plt.tight_layout() # looks strange otherwise
plt.show()



# now, we're looking at antibiotics exposure. the original AbxPreCollection 
# listed different types of antibiotics and also had the option of no_abx. 
# so now we turn this into a binary column, whose values are true or false 
# based on if we've had abx
df['AbxPreCollection_binary'] = df['AbxPreCollection'] != 'no_abx' # true if abx esposure, false if not

# these are the nodes we're interested in
corr_nodes = ['38', '51', '57', '18', '65', '76', '88']
corr_nodes.sort()  # so they display better

# rearrange dataframe
df_melted = pd.melt(df, id_vars='AbxPreCollection_binary', value_vars=corr_nodes, var_name='vals', value_name='Value')

# to get all of the text to be visible
plt.rcParams.update({
    'axes.titlesize': 18,
    'axes.labelsize': 16,
    'xtick.labelsize': 14,
    'ytick.labelsize': 14,
    'legend.fontsize': 14,
    'legend.title_fontsize': 15
})

# make the stripplot
plt.figure(figsize=(10, 6))
ax = sns.stripplot(x='vals',
                   y='Value',
                   hue='AbxPreCollection_binary', # color by whether or not abx
                   data=df_melted,
                   jitter=True, alpha=0.3, size=6)

plt.title('Node Values colored by Abx (Binary)')
plt.xlabel('Node')
plt.ylabel('Value')
plt.legend(title='Abx')

# annotate each node with the p-value
for i, col in enumerate(corr_nodes):
    # run the t-test
    group_true = df[df["AbxPreCollection_binary"] == True][col]
    group_false = df[df["AbxPreCollection_binary"] == False][col]
    t_stat, p_value = ttest_ind(group_true, group_false, nan_policy='omit')

    # figure out where to position the annotation
    max_val = df[col].max()
    y_pos = max_val + 0.05 * (df_melted['Value'].max() - df_melted['Value'].min())

    # format p-value
    p_text = f"p = {p_value:.2e}"
    # plot it
    plt.text(i, y_pos, p_text, ha='center')

plt.tight_layout()
plt.show()
