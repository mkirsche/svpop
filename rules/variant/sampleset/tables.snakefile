"""
Tables of sampleset variants.
"""

SVTYPE_ORDER = ['INS', 'DEL', 'INV', 'SNV']


# variant_tables_caller_xlsx
#
# TSV to Excel.
rule variant_sampleset_table_unmerged_variants_xlsx:
    input:
        tab='results/variant/sampleset/{sampleset}/tables/{samplelist}/{svset}/{filter}/pre_merge/{vartype}_{svtype}.tab'
    output:
        xlsx='results/variant/sampleset/{sampleset}/tables/{samplelist}/{svset}/{filter}/pre_merge/{vartype}_{svtype}.xlsx'
    run:

        pd.read_csv(
            input.tab, sep='\t'
        ).to_excel(
            output.xlsx, index=False
        )

rule variant_sampleset_table_unmerged_variants:
    input:
        tab=lambda wildcards: analib.sampleset.get_sample_set_input(
            wildcards.sampleset,
            wildcards.samplelist,
            'results/variant/{sourcetype}/{sourcename}/tables/{sample}/{svset}/{filter}/variant_summary/{vartype}_{svtype}.tab',
            config,
            wildcards
        )
    output:
          tab='results/variant/sampleset/{sampleset}/tables/{samplelist}/{svset}/{filter}/pre_merge/{vartype}_{svtype}.tab'
    run:

        sampleset_entry = analib.sampleset.get_config_entry(wildcards.sampleset, wildcards.samplelist, config)

        # Read
        df = pd.concat(
            [pd.read_csv(file_name, sep='\t') for file_name in input.tab],
            axis=0
        )

        if 'SOURCETYPE' in df:
            del(df['SOURCETYPE'])

        if 'SOURCENAME' in df:
            del(df['SOURCENAME'])

        # Get a list of columns and SVTYPE values keeping them ordered.
        col_names = [val for val in df.columns if val not in {'SAMPLE', 'SVTYPE'}]

        svtype_list = list()

        for val in df['SVTYPE']:
            if val not in svtype_list:
                svtype_list.append(val)

        svtype_list_head = [val for val in SVTYPE_ORDER if val in svtype_list]
        svtype_list_tail = [val for val in svtype_list if val not in svtype_list_head]

        svtype_list = svtype_list_head + svtype_list_tail

        # Reshape: Multi-index of values in svtype_list and col_names
        df = df.pivot(index='SAMPLE', columns='SVTYPE')

        df = df.swaplevel(0, 1, axis=1)

        df = df.reindex([(svtype, col) for svtype in svtype_list for col in col_names], axis=1)

        # Order samples
        df = df.reindex(sampleset_entry['samples'])

        # Save
        df.to_csv(output.tab, sep='\t')
