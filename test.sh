cp constraints/* models/* rates/* medium/* .


. ./adynDFBA.R -m iJV1220 -c wt_glc+aa -n NMMP+glc -s 30 -V 5 -a
. ./adynDFBA.R -m iJV1220 -c plasmid_glc+aa -n NMMP+glc -s 30 -V 5 -a
. ./adynDFBA.R -m iJV1220 -c tnf_glc+aa -n NMMP+glc -s 30 -V 5 -a

. ./adynDFBA.R -m iJV1220 -c aml_mnl+aa -n NMMP+mnl -r rates.aml -s 30 -V 5 -a
. ./adynDFBA.R -m iJV1220 -c dag_mnl+aa -n NMMP+mnl -r rates.dag -s 70 -V 5 -a

rm *.tab *.tsv *.dat
