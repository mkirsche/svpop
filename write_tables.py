import pandas as pd
import sys
import numpy as np
from analib.variant import get_variant_id
from analib.variant import order_variant_columns

fn = sys.argv[1]

df = pd.read_csv(fn, sep='\t', header=0)

# Remove uninformative columns
del(df['REF'], df['ALT'], df['QUAL'])

# Normalize fields
df['SVLEN'] = np.abs(df['SVLEN'])

df['SVTYPE_ORG'] = df['SVTYPE']
df['SVTYPE'] = df['SVTYPE'].apply(lambda val: val.replace('/', ''))

df['END'] = df.apply(lambda row: row['POS'] + (row['SVLEN'] if row['SVTYPE'] != 'INS' else 1), axis=1)

# Set ID and arrange columns
df['ID2'] = get_variant_id(df)

df = order_variant_columns(df)

# Sort
# Sniffles sometimes writes duplicate entries for variant calls. Sorting by GT (descending so 1/1 comes first)
# DV (descending so higher values of supporting reads come first) and DR (ascending so lower values reference
# reads come first) allows duplicates to be removed so that keeping the first variant in the table gives the
# strongest version of the variant.
df.sort_values(
    ['#CHROM', 'POS', 'SVTYPE', 'SVLEN', 'GT', 'DV', 'DR'],
    ascending=(True, True, True, True, False, False, True),
    inplace=True
)

# Write full caller annotations
#df.to_csv(output.bed_anno, sep='\t', index=False)

# Filter by VCF column
df = df.loc[df['FILTER'] == 'PASS']

# Drop duplicates (Sniffles sometimes writes duplicate variant calls).
df.drop_duplicates(['ID2'], inplace=True)

# Remove annotation columns (recorded in the annotation table already written)
for rm_field in (
    'PRECISE', 'IMPRECISE',
    'Kurtosis_quant_start', 'Kurtosis_quant_stop', 'STD_quant_start', 'STD_quant_stop',
    'MAPQ', 'RE', 'REF_strand', 'STRANDS', 'SUPTYPE', 'AF', 'SVMETHOD', 'ZMW', 'SAMPLE',
    'FILTER', 'SVTYPE_ORG'
):
    del(df[rm_field])


# Rename fields
prepend_caller_set = {'GT', 'DV', 'DR', 'CHR2'}

df.columns = [('SNIF_' + col if col in prepend_caller_set else col) for col in df.columns]

# Get SV tables
df_ins = df.loc[df['SVTYPE'] == 'INS'].copy()
df_del = df.loc[df['SVTYPE'] == 'DEL'].copy()
df_inv = df.loc[df['SVTYPE'] == 'INV'].copy()
df_dup = df.loc[(df['SVTYPE'] == 'DUP') | (df['SVTYPE'] == 'INVDUP')].copy()
df_bnd = df.loc[df['SVTYPE'] == 'BND'].copy()

# Annotate inverted duplications
df_dup['SNIF_INVDUP'] = df_dup['SVTYPE'] == 'INVDUP'
df_dup['SVTYPE'] = 'DUP'

# CHR2 only has information for BND and INS
del(df_del['SNIF_CHR2'])
del(df_inv['SNIF_CHR2'])
del(df_dup['SNIF_CHR2'])

df_ins.loc[df_ins['SVLEN'] >= 50].to_csv(fn + ".sv_ins", sep='\t', index=False)
df_del.loc[df_del['SVLEN'] >= 50].to_csv(fn + ".sv_del", sep='\t', index=False)
df_inv.loc[df_inv['SVLEN'] >= 50].to_csv(fn + ".sv_inv", sep='\t', index=False)
df_dup.loc[df_dup['SVLEN'] >= 50].to_csv(fn + ".sv_dup", sep='\t', index=False)
df_bnd.to_csv(fn + ".sv_bnd", sep='\t', index=False)
