
from hail import *
hc = HailContext()

sex = sys.argv[1]             # {'both_sexes', 'female', 'male'}
pipeline_type = sys.argv[2]   # {'phesant', 'icd10', 'finngen'}
try:
    pipeline_number = sys.argv[3] # {0, 1, 2, ...}
except:
    pipeline_number = ''

if sex not in set(['both_sexes', 'female', 'male']):
    raise ValueError('Invalid sex value "{}" - must be one of {"both_sexes", "female", "male"}.'.format(sex))
if pipeline_type not in set(['phesant', 'icd10', 'finngen']):
    raise ValueError('Invalid pipeline type "{}" - must be one of {"phesant", "icd10", "finngen"}.'.format(pipeline_type))

print '#####################'
print '## Starting... ######'
print '## Sex: {}'.format(sex)
print '## Pipeline type: {}'.format(pipeline_type)
print '## Pipeline number: {}'.format(pipeline_number)
print '#####################'

def group_phenotypes(kt, block_type):
    codes = [x for x in kt.columns if x != 's']
    if block_type == 'single':
        phenotype_groups = [['single_block', codes]]
    elif block_type == 'individual':
        phenotype_groups = [[c, [c]] for c in codes]
    elif block_type == 'mix':
        phenotype_groups = []
        for code in codes:
            prefix = code.split('_')[0]
            try:
                idx = [x[0] for x in phenotype_groups].index(prefix)
            except ValueError:
                phenotype_groups.append([prefix, [code]])
            else:
                phenotype_groups[idx][1].append(code)
    return phenotype_groups

pipeline_path = 'gs://ukb31063-mega-gwas/phenotype-pipelines/{0}/ukb31063.{1}.{0}.pipeline.{2}.kt'.format(pipeline_type, sex, pipeline_number)
results_auto_path = 'gs://ukb31063-mega-gwas/dominance-results-tables/{0}/ukb31063.{1}.{0}.autosomes.dominance.pipeline.{2}.kt'.format(pipeline_type, sex, pipeline_number)
results_chrX_path = 'gs://ukb31063-mega-gwas/dominance-results-tables/{0}/ukb31063.{1}.{0}.chrX.dominance.pipeline.{2}.kt'.format(pipeline_type, sex, pipeline_number)

kt_pipeline = hc.read_table(pipeline_path)
kt_results_auto = hc.read_table(results_auto_path)
kt_results_chrX = hc.read_table(results_chrX_path)
kt_results = KeyTable.union(kt_results_auto, kt_results_chrX)
kt_results = kt_results.annotate('variant = Variant(v.contig.replace("^0", ""), v.start, v.ref, v.alt())')
kt_results = (kt_results.key_by('variant')
                        .order_by('variant')
                        .filter('isDefined(va.AF)')
                        .drop('v')).cache()

if pipeline_type == 'phesant':
    kt_phenosummary = hc.import_table('gs://ukb31063-mega-gwas/phenotype-summaries/phesant/ukb31063.{0}.phenosummary.*.tsv'.format(sex)).cache()
    phenotype_groups = group_phenotypes(kt_pipeline, block_type='mix')
    kt_pipeline = hc.read_table('gs://ukb31063-mega-gwas/phenotype-pipelines/phesant/ukb31063.{0}.phesant.pipeline.{1:}.kt'.format(sex, pipeline_number))
elif (pipeline_type == 'icd10') or (pipeline_type == 'finngen'):
    kt_phenosummary = hc.import_table('gs://ukb31063-mega-gwas/phenotype-summaries/{0}/ukb31063.{1}.{0}.phenosummary.pipeline.{2:}.tsv'.format(pipeline_type, sex, pipeline_number)).cache()
    phenotype_groups = group_phenotypes(kt_pipeline, block_type='single')

count = 1
for i, group in enumerate(phenotype_groups):
    for j, code in enumerate(group[1]):
        print 'Exporting {0} ({1:})...'.format(code, count)
        kt_export = kt_results.annotate(['n_complete_samples = va.results[{0:}].nCompleteSamples'.format(i),
                                         'AC = va.results[{0:}].AC'.format(i),
                                         'ytx = va.results[{0:}].ytx[{1:}]'.format(i, j),
                                         'beta = va.results[{0:}].beta[{1:}]'.format(i, j),
                                         'se = va.results[{0:}].se[{1:}]'.format(i, j),
                                         'tstat = va.results[{0:}].tstat[{1:}]'.format(i, j),
                                         'pval = va.results[{0:}].pval[{1:}]'.format(i, j)])
        kt_export = kt_export.annotate('AF = AC/(2.0 * n_complete_samples)')
        kt_export = kt_export.annotate(['minor_allele = if (AF <= 0.5) variant.alt() else variant.ref',
                                        'minor_AF = if (AF <= 0.5) AF else 1.0 - AF'])
        if pipeline_type == 'phesant':
            try:
                n_cases = int(kt_phenosummary.query('`N.cases`.filter(x => FieldID == "{}").collect()'.format(code))[0])
            except TypeError:
                value_counter = kt_pipeline.query('`{}`.map(x => str(x)).counter()'.format(code))
                if len(value_counter) <= 4:
                    min_category_count = min([int(x) for x in value_counter.values()])
                    kt_export = kt_export.annotate('expected_min_category_minor_AC = 2.0 * minor_AF * {:}.toInt()'.format(min_category_count))
                    kt_export = kt_export.annotate('low_confidence_variant = (expected_min_category_minor_AC < 25) || (minor_AF < 0.001)')
                else:
                    kt_export = kt_export.annotate('low_confidence_variant = minor_AF < 0.001')
            else:
                kt_export = kt_export.annotate('expected_case_minor_AC = 2.0 * minor_AF * {:}.toInt()'.format(n_cases))
                kt_export = kt_export.annotate('low_confidence_variant = (expected_case_minor_AC < 25) || (minor_AF < 0.001)')

        elif pipeline_type == 'icd10' or pipeline_type == 'finngen':
            n_cases = kt_phenosummary.query('n_cases.filter(x => code == "{}").collect()'.format(code))[0]
            kt_export = kt_export.annotate('expected_case_minor_AC = 2.0 * minor_AF * {:}.toInt()'.format(n_cases))
            kt_export = kt_export.annotate('low_confidence_variant = (expected_case_minor_AC < 25) || (minor_AF < 0.001)')
        
        try:
            kt_export = kt_export.select(['variant',
                                          'minor_allele',
                                          'minor_AF',
                                          'expected_case_minor_AC',
                                          'low_confidence_variant',
                                          'n_complete_samples',
                                          'AC',
                                          'ytx',
                                          'beta',
                                          'se',
                                          'tstat',
                                          'pval'])
        except:
            kt_export = kt_export.select(['variant',
                                          'minor_allele',
                                          'minor_AF',
                                          'low_confidence_variant',
                                          'n_complete_samples',
                                          'AC',
                                          'ytx',
                                          'beta',
                                          'se',
                                          'tstat',
                                          'pval'])

        kt_export.export('gs://ukb31063-mega-gwas/dominance-results-tsvs/{0}.dominance.gwas.imputed_v3.{1}.tsv.bgz'.format(code, sex))
        count += 1

print '#####################'
print '## COMPLETED ######'
print '## Sex: {}'.format(sex)
print '## Pipeline type: {}'.format(pipeline_type)
print '## Pipeline number: {}'.format(pipeline_number)
print '#####################'
