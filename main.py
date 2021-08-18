import argparse
from phaser import Phaser
from family import Family

#!/usr/bin/env python3
"""
Phase variants in a VCF with the WhatsHap algorithm

Read a VCF and one or more files with phase information (BAM/CRAM or VCF phased
blocks) and phase the variants. The phased VCF is written to standard output.
"""
import logging
import sys
import platform
from collections import defaultdict
from contextlib import ExitStack
from typing import Optional, List, TextIO
from vcf import VcfError
from graph import ComponentFinder
from pedigree import (
    PedReader,
    mendelian_conflict,
)
import ped
import ped_utils
from ped_utils import Trio
from timer import StageTimer
from utils import plural_s, warn_once, NumericSampleIds
from phasevcf import CommandLineError, log_memory_usage, PhasedInputReader

logger = logging.getLogger(__name__)


def setup_pedigree(ped_path, samples):
    """
    Read in PED file to set up list of relationships.

    Return a pair (trios, pedigree_samples), where trios is a list of Trio
    objects and pedigree_samples is the set of all samples that are mentioned
    in the PED file (as individual, mother or father).

    ped_path -- path to PED file
    samples -- samples to be phased
    """
    trios = []
    pedigree_samples = set()
    for trio in PedReader(ped_path):
        if trio.child is None or trio.mother is None or trio.father is None:
            warn_once(
                logger,
                "Relationship %s/%s/%s ignored because at least one of the individuals is unknown.",
                trio.child,
                trio.mother,
                trio.father,
            )
            continue
        # if at least one individual is not in samples, skip trio
        if (
            (trio.mother not in samples)
            or (trio.father not in samples)
            or (trio.child not in samples)
        ):
            # happens in case --ped and --samples are used
            warn_once(
                logger,
                "Relationship %s/%s/%s ignored because at least one of the "
                "individuals was not given by --samples.",
                trio.child,
                trio.mother,
                trio.father,
            )
            continue

        trios.append(trio)
        pedigree_samples.add(trio.child)
        pedigree_samples.add(trio.father)
        pedigree_samples.add(trio.mother)

    return trios, pedigree_samples


def phase_s1_with_s2(s1, s2):
    pass


def run_pedhap(
    variant_file: str,
    output: TextIO = sys.stdout,
    chromosomes: Optional[List[str]] = None,
    indels: bool = True,
    mapping_quality: int = 20,
    max_coverage: int = 15,
    ped: Optional[str] = None,
    recombrate: float = 1.26,
    genmap: Optional[str] = None,
    genetic_haplotyping: bool = True,
    default_gq: int = 30,
    use_ped_samples: bool = True,
    tag: str = "PS",
):
    timers = StageTimer()
    logger.info(
        f"This is pedhap , running under Python {platform.python_version()}")
    numeric_sample_ids = NumericSampleIds()
    command_line: Optional[str]
    command_line = None
    with ExitStack() as stack:
        try:
            # vcf_writer = stack.enter_context(
            #     PhasedVcfWriter(
            #         command_line=command_line,
            #         in_path=variant_file,
            #         out_file=output,
            #         tag=tag,
            #         indels=indels,
            #     )
            # )
            pass
        except (OSError, VcfError) as e:
            raise CommandLineError(e)

        input_vcf = PhasedInputReader(
            variant_file, indels=indels,)
        show_phase_vcfs = input_vcf.has_vcfs

        # if --use-ped-samples is set, use only samples from PED file
        if ped and use_ped_samples:
            samples = PedReader(ped).samples()

        raise_if_any_sample_not_in_vcf(input_vcf.vcf_reader, samples)

        families, family_trios = setup_families(samples, ped, max_coverage)
        del samples
        for trios in family_trios.values():
            for trio in trios:
                # Ensure that all mentioned individuals have a numeric id
                _ = numeric_sample_ids[trio.child]

        with timers("parse_phasing_vcfs"):
            # TODO should this be done in PhasedInputReader.__init__?
            input_vcf.read_vcfs()

        for chromosome, variant_table in input_vcf.chromVaritables.items():
            if (not chromosomes) or (chromosome in chromosomes):
                logger.info("======== Working on chromosome %r", chromosome)
            else:
                logger.info(
                    "Leaving chromosome %r unchanged (present in VCF but not requested by option --chromosome)",
                    chromosome,
                )
                # vcf_writer.write(chromosome)

            # Iterate over all families to process, i.e. a separate DP table is created
            # for each family.
            # TODO: Can the body of this loop be factored out into a phase_family function?
            for representative_sample, family in sorted(families.items()):
                logger.info("---- Processing individual %s",
                            representative_sample)
                max_coverage_per_sample = max(1, max_coverage // len(family))
                logger.info("Using maximum coverage per sample of %dX",
                            max_coverage_per_sample)
                trios = family_trios[representative_sample]
                len(trios) > 0
                s1 = trio.child
                s2 = trio.father
                s3 = trio.mother
                print(s1)
                print(s2)
                print("xxx")
                variant_table.phase_single(s1, s2, side=0)
                variant_table.phase_single(s1, s3, side=1)
                variant_table.write(s1, output)
                break
                # vcf_writer.write(chromosome)
                # for sample in family:
                #     trios
            logger.debug("Chromosome %r finished", chromosome)

    log_time_and_memory_usage(timers, show_phase_vcfs=show_phase_vcfs)


def raise_if_any_sample_not_in_vcf(vcf_reader, samples):
    vcf_sample_set = set(vcf_reader.samples)
    for sample in samples:
        if sample not in vcf_sample_set:
            raise CommandLineError(
                "Sample {!r} requested on command-line not found in VCF".format(
                    sample)
            )


def setup_families(samples, ped, max_coverage):
    """
    Return families, family_trios pair.

    families maps a family representative to a list of family members

    family_trios maps a family representative to a list of trios in this family
    """

    # list of all trios across all families
    all_trios = dict()

    # Keep track of connected components (aka families) in the pedigree
    family_finder = ComponentFinder(samples)

    if ped:
        all_trios, pedigree_samples = setup_pedigree(ped, samples)

        for trio in all_trios:
            family_finder.merge(trio.father, trio.child)
            family_finder.merge(trio.mother, trio.child)

    # map family representatives to lists of family members
    families = defaultdict(list)
    for sample in samples:
        families[family_finder.find(sample)].append(sample)

    # map family representatives to lists of trios for this family
    family_trios = defaultdict(list)
    for trio in all_trios:
        family_trios[family_finder.find(trio.child)].append(trio)
    logger.info(
        "Working on %d%s samples from %d famil%s",
        len(samples),
        plural_s(len(samples)),
        len(families),
        "y" if len(families) == 1 else "ies",
    )

    largest_trio_count = max([0] + [len(trio_list)
                             for trio_list in family_trios.values()])
    if max_coverage + 2 * largest_trio_count > 23:
        logger.warning(
            "The maximum coverage is too high! "
            "WhatsHap may take a long time to finish and require a huge amount of memory."
        )
    return families, family_trios


def log_time_and_memory_usage(timers, show_phase_vcfs):
    total_time = timers.total()
    logger.info("\n== SUMMARY ==")
    log_memory_usage()
    # fmt: off
    # logger.info("Time spent reading BAM/CRAM:                 %6.1f s", timers.elapsed("read_bam"))
    # logger.info("Time spent parsing VCF:                      %6.1f s", timers.elapsed("parse_vcf"))
    # if show_phase_vcfs:
    #     logger.info("Time spent parsing input phasings from VCFs: %6.1f s", timers.elapsed("parse_phasing_vcfs"))
    # logger.info("Time spent selecting reads:                  %6.1f s", timers.elapsed("select"))
    # logger.info("Time spent phasing:                          %6.1f s", timers.elapsed("phase"))
    # logger.info("Time spent writing VCF:                      %6.1f s", timers.elapsed("write_vcf"))
    # logger.info("Time spent finding components:               %6.1f s", timers.elapsed("components"))
    # logger.info("Time spent on rest:                          %6.1f s", total_time - timers.sum())
    # logger.info("Total elapsed time:                          %6.1f s", total_time)
    # fmt: on


def find_mendelian_conflicts(trios, variant_table):
    mendelian_conflicts = set()
    for trio in trios:
        genotypes_mother = variant_table.genotypes_of(trio.mother)
        genotypes_father = variant_table.genotypes_of(trio.father)
        genotypes_child = variant_table.genotypes_of(trio.child)

        for index, (gt_mother, gt_father, gt_child) in enumerate(
            zip(genotypes_mother, genotypes_father, genotypes_child)
        ):
            if (not gt_mother.is_none()) and (not gt_father.is_none()) and (not gt_child.is_none()):
                if mendelian_conflict(gt_mother, gt_father, gt_child):
                    mendelian_conflicts.add(index)
    return mendelian_conflicts

# def phasing_trio_child(phaser: Phaser, trio: Trio, chromo: str):
#     phaser.phasing_trio_child(chromo, trio)


def up_to_down(all_trios: List[Trio], phaser: Phaser, chromo):
    # find top level trio and phsing child
    top_level_trios = ped_utils.get_top_level_trios(all_trios)
    # t = top_level_trios[0]
    # phaser.phasing_trio_child(t)
    for t in top_level_trios:
        phaser.phasing_trio_child(t, chromo)
    next_level_trios = ped_utils.get_next_level_trios(
        all_trios, top_level_trios)
    while len(next_level_trios) != 0:
        for t in next_level_trios:
            phaser.phasing_trio_child(t, chromo)
        next_level_trios = ped_utils.get_next_level_trios(
            all_trios, next_level_trios)


def down_to_up(all_trios: List[Trio], phaser: Phaser, chromo):
    # find bottom level trio and phsing child
    bottom_level_trios = ped_utils.get_bottom_level_trios(all_trios)
    # t = top_level_trios[0]
    # phaser.phasing_trio_child(t)
    # phaser.write_phased_result(t.child.id, "/home/caronkey/Documents/cityu/pedhap/test/test2.out")
    for t in bottom_level_trios:
        phaser.phasing_trio_parent(t, chromo)
    prev_level_trios = ped_utils.get_prev_level_trios(
        all_trios, bottom_level_trios)
    while len(prev_level_trios) != 0:
        for t in prev_level_trios:
            phaser.phasing_trio_parent(t, chromo)
        prev_level_trios = ped_utils.get_prev_level_trios(
            all_trios, prev_level_trios)


def iter_phase(all_trios: List[Trio], phaser: Phaser):
    up_to_down(all_trios, phaser)
    down_to_up(all_trios, phaser)


def main():
    parser = argparse.ArgumentParser("trio phase")
    parser.add_argument(
        '-v', help='merged VCF file', required=False, dest='vcf_file')
    parser.add_argument(
        '-p', help='pedigree file', required=False, dest='ped_file')
    parser.add_argument(
        '-o', help='out phased vcf file', required=False, dest='out_file')
    parser.add_argument(
        '--max_round', help='max phasing iter times, if not given, decided by program', 
        required=False, default=0, type=int,
        dest='max_round')
    args = parser.parse_args()

    phaser = Phaser(vcf_file=args.vcf_file, out_file=args.out_file, max_round = int(args.max_round))
    families = ped.open_ped(args.ped_file)
    for f in families:
        all_trios = ped_utils.get_trios(f)
        for chromo in phaser.chromos:
            while phaser.check_phasing_state(chromo):
                up_to_down(all_trios, phaser, chromo)
                down_to_up(all_trios, phaser, chromo)
        # iter_phase(all_trios, phaser, args.max_round)
    # phaser.write()
    phaser.write_simple("s0210-1_FDHG190451805-1a")
        # phaser.write_phased_result("s0210-1_FDHG190451805-1a", "/home/caronkey/Documents/cityu/pedhap/test/test2.out")
    # run_pedhap(variant_file=args.vcf_file,
    #            ped=args.ped_fle, output=args.out_file)


if __name__ == "__main__":
    main()
