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
    bush = subparsers.add_parser(
        "bush-gui", help="launch the bush clearance / rotor-mouth animation (coupled orbit)"
    )
    bush.add_argument("--gif", type=str, default=None, help="render a GIF to this path and exit")
    bush.add_argument("--cache", type=str, default=None, help="orbit .npz cache (load or write)")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """Run the command-line interface."""

    parser = build_parser()
    arguments = parser.parse_args(argv)
    if arguments.command == "gui":
        from mochi.gui import main as run_gui

        return run_gui()
    if arguments.command == "bush-gui":
        from mochi.bush_gui import main as run_bush_gui

        forward: list[str] = []
        if arguments.gif:
            forward += ["--gif", arguments.gif]
        if arguments.cache:
            forward += ["--cache", arguments.cache]
        return run_bush_gui(forward)
    parser.print_help()
    return 0
