from flytekit import workflow
from run.tests.test_index_cram import test_index_cram_wf
from run.tests.test_split_cram import test_split_cram_wf
from run.tests.test_golem_call_variants import test_golem_call_variants_wf

@workflow
def registered():
    test_index_cram_wf()
    test_split_cram_wf()
    test_golem_call_variants_wf()