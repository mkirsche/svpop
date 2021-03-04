from analib.svmerge import merge_variants
import pandas as pd
import sys

filenames=[]
filelist=sys.argv[1]
with open(filelist) as f:
	for line in f:
		if len(line) > 0:
			filenames.append(line.strip())
print(filenames)

samples = sys.argv[2].split(",")
print(samples)

outprefix = sys.argv[3]
print(outprefix)

fulldf = pd.DataFrame()
types = ['ins', 'del', 'inv', 'dup', 'bnd']
for type in types:
	newlist = []
	for filename in filenames:
		newlist.append(filename + ".sv_" + type)
	df = merge_variants(newlist, samples, 'nr:szro=50:offset=200', threads=1)

	# Bylen to byref
	if df.shape[0] > 0:
		df['END'] = df.apply(lambda row: (row['POS'] + 1) if row['SVTYPE'] == 'INS' else row['END'], axis=1)

	# Write
	df.to_csv(outprefix + 'merged_' + type + '.bed', sep='\t', index=False)
	fulldf = pd.concat([fulldf, df])
	
fulldf.to_csv(outprefix + 'merged_all.bed', sep='\t', index=False)
