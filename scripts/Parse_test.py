import argparse

def argument_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--metag",
                        action="store_true",
                        help="Use to activate function relevant to metagenomics pipeline")
    return parser


args = argument_parser().parse_args()

if args.metag:
    print("Positional written")
else:
    print("Positional not there")
