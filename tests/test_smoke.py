from pytest import CaptureFixture

from mochi import __version__
from mochi.cli import build_parser, main


def test_package_has_version() -> None:
    assert __version__ == "0.1.0"


def test_cli_without_command_shows_help(capsys: CaptureFixture[str]) -> None:
    assert main([]) == 0
    assert "usage: mochi" in capsys.readouterr().out


def test_gui_command_is_registered_without_launching_window() -> None:
    arguments = build_parser().parse_args(["gui"])
    assert arguments.command == "gui"
