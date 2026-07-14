"""Command-line entry point for mochi."""

from argparse import ArgumentParser
from collections.abc import Sequence

from mochi import __version__


def build_parser() -> ArgumentParser:
    """Build the command-line parser."""

    parser = ArgumentParser(
        prog="mochi",
        description="Rotary-compressor simulation utilities.",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    subparsers = parser.add_subparsers(dest="command")
    subparsers.add_parser("gui", help="launch the supporting motion visualization")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """Run the command-line interface."""

    parser = build_parser()
    arguments = parser.parse_args(argv)
    if arguments.command == "gui":
        from mochi.gui import main as run_gui

        return run_gui()
    parser.print_help()
    return 0
