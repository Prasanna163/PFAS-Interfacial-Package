import sys

from pfas_interface_cli.workflow import build_parser, main, run_workflow


__all__ = ["main", "build_parser", "run_workflow"]


if __name__ == "__main__":
    main(sys.argv[1:])
