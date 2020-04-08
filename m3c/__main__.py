import argparse
import functools
import logging
import sys

from m3c import prefill
import metab_admin as server
import metab_import as generate
import pubfetch


PROGRAM: str = "m3c"

logger = logging.getLogger(PROGRAM)


def init_logger(log: str, verbose: bool):
    # Create the file logger.
    logfile = logging.FileHandler(f"{PROGRAM}.log", delay=True)
    logfile.setFormatter(logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    ))

    # Create a console logger.
    console = logging.StreamHandler(sys.stderr)
    if not verbose:
        console.setLevel(logging.INFO)
    # Use a terser format for a user interacting with our program.
    console.setFormatter(logging.Formatter(
        f"%(relativeCreated)15f %(levelname)+{len('CRITICAL')}s  %(message)s"
    ))

    logging.basicConfig(level=logging.DEBUG, handlers=[console, logfile])


def main():
    args = parse_args(sys.argv[1:])
    init_logger(args.log, args.verbose)

    logger.debug(f"{PROGRAM} {args.cmd} started")

    if args.cmd == "prefill":
        prefill.prefill(args.config)
    elif args.cmd == "serve":
        server.serve(args.config)
    elif args.cmd == "generate":
        generate.generate(args.config, args.add_devs, args.diff)
    elif args.cmd == "pubfetch":
        pubfetch.pubfetch(args.config, args.authorships, args.delay, args.max)

    logger.debug(f"{PROGRAM} ended")


def parse_args(args):
    parser = argparse.ArgumentParser(prog=PROGRAM)
    nat = functools.partial(natural, parser)
    parser.add_argument("-l", "--log", default=f"{PROGRAM}.log",
                        help="path to the log file")
    parser.add_argument("-v", "--verbose", action="store_true", default=False,
                        help="increase verbosity of console log")
    subparsers = parser.add_subparsers(
        title="Sub-commands", description="Valid subcommands", dest="cmd"
    )
    subparsers.add_parser(
        "prefill",
        help="downloads and processes data from the Metabolomics Workbench",
    )
    generate = subparsers.add_parser(
        "generate",
        help="generates N-Triples for import into the People Portal"
    )
    generate.add_argument(
        "-a", "--add-devs", action="store_true", default=True, dest="add_devs",
        help="add a Person record for new tools developers"
    )
    generate.add_argument("-x", "--diff", help="path for differential update")
    pubfetchcmd = subparsers.add_parser(
        "pubfetch", help="downloads PubMed publication data"
    )
    pubfetchcmd.add_argument(
        "--authorships", action="store_true", default=False,
        help="only update authorships; do not download publications"
    )
    pubfetchcmd.add_argument(
        "--delay", type=nat, default=0,
        help="number of seconds to wait between PubMed requests"
    )
    pubfetchcmd.add_argument(
        "--max", type=nat, default=None,
        help="maximum number of authorship searches to perform"
    )
    subparsers.add_parser(
        "serve",
        help="starts an HTTP server for the Admin Forms"
    )
    parser.add_argument("config", help="path to the YAML configuration file")
    parsed = parser.parse_args(args)
    return parsed


def natural(parser: argparse.ArgumentParser, value) -> int:
    i = int(value)
    if i < 0:
        parser.exit(2, f"value must be non-negative: got {value}\n")
    return i


if __name__ == "__main__":
    main()
