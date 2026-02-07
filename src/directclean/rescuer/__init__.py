"""
Internal adapter rescuer for Direct-cDNA reads.

Detects reads where two (or more) cDNA molecules were ligated
together via internal TSO/RTP adapters, then chops them into
independent sub-reads that can be aligned separately.

Public API::

    from directclean.rescuer import ReadChopper, AdapterConfig

    chopper = ReadChopper(
        input_fastq="restranded.fastq",
        output_fastq="rescued.fastq",
        config=AdapterConfig(max_edit_distance=3),
    )
    report = chopper.run()
"""

from directclean.rescuer.adaptor_seq import (
    AdapterConfig,
    TSO_SEQUENCE,
    RTP_SEQUENCE,
    RTP_RC_SEQUENCE,
)
from directclean.rescuer.adapter_finder import (
    AdapterFinder,
    AdapterHit,
    InternalJunction,
    FinderResult,
)
from directclean.rescuer.chopper import (
    ReadChopper,
    RescueReport,
)

__all__ = [
    "AdapterConfig",
    "TSO_SEQUENCE",
    "RTP_SEQUENCE",
    "RTP_RC_SEQUENCE",
    "AdapterFinder",
    "AdapterHit",
    "InternalJunction",
    "FinderResult",
    "ReadChopper",
    "RescueReport",
]