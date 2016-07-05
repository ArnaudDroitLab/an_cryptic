export RAP_ID="eav-760-aa"

mkdir output
$MUGQIC_PIPELINES_HOME/pipelines/rnaseq/rnaseq.py -s '1-22' -l debug \
    -r input/readset.txt \
    -d input/design.txt \
    -o /gs/project/eav-760-aa/RNA-CKII-Spt6/output-new \
    --config ./rnaseq.base.ini \
        $MUGQIC_PIPELINES_HOME/pipelines/rnaseq/rnaseq.guillimin.ini \
        $MUGQIC_PIPELINES_HOME/resources/genomes/config/Saccharomyces_cerevisiae.R64-1-1.ini \
        ./this_run.ini