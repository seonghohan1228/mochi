"""Allow ``python -m mochi`` to run the command-line interface."""

from mochi.cli import main

raise SystemExit(main())
