"""Command-line entry point for mochi."""

from argparse import ArgumentParser
from collections.abc import Sequence

from mochi import __version__


def build_parser() -> ArgumentParser:
    """Build the command-line parser."""

    parser = ArgumentParser(
        prog="mochi",
        description="Numerical solver for scroll-compressor bushing motion.",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """Run the command-line interface."""

    parser = build_parser()
    parser.parse_args(argv)
    parser.print_help()
    return 0
