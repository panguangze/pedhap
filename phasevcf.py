import sys
import resource
import logging
from vcf import VcfReader

logger = logging.getLogger(__name__)


class CommandLineError(Exception):
    """An anticipated command-line error occurred. This ends up as a user-visible error message"""

class PhasedInputReader:
    def __init__(
        self,
        vcf_path,
        indels,
        **kwargs,  # passed to ReadSetReader constructor
    ):
        self._vcf_path = vcf_path
        # TODO exit stack!
        vcf_reader = VcfReader(path = self._vcf_path, indels=indels, phases=True)
        self.vcf_reader = vcf_reader
        self.chromVaritables = {}

    def __enter__(self):
        return self

    def __exit__(self, *args):
        if self._fasta is not None:
            self._fasta.close()

    @property
    def has_vcfs(self):
        return bool(self._vcf_path)

    def read_vcfs(self):
        # Read phase information provided as VCF files, if provided.
        # TODO: do this chromosome- and/or sample-wise on demand to save memory.
        logger.info("Reading phased blocks from %r", self.vcf_reader.path)
        for variant_table in self.vcf_reader:
            self.chromVaritables[variant_table.chromosome] = variant_table

def log_memory_usage(include_children=False):
    if sys.platform == "linux":
        if include_children:
            memory_kb = (
                resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
                + resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
            )
        else:
            memory_kb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        logger.info("Maximum memory usage: %.3f GB", memory_kb / 1e6)
