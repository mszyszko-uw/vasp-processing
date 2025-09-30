import argparse

avalable_steps = ('dry', 'scf', 'geo', 'bs', 'dos', 'eps', 'epsfcn', 'meff', 'selrul', 'gband', 'gexc', 'test_1', 'test_2')

def parse_part_values(values):
    if len(values) == 1 and values[0].lower() == 'all':
        return avalable_steps

    if any(val.lower() == 'all' for val in values):
        raise argparse.ArgumentTypeError("Use only 'all' value")
    result = set()
    steps_to_run = list( dict.fromkeys( values.lower().split(',') )) 

    if not all(x in avalable_steps for x in steps_to_run):
        raise ValueError('Unknow step')
    return steps_to_run

def print_avalable_steps():
    print(avalable_steps)

def read_cli_params():
#    parser = argparse.ArgumentParser(description="Usage: script.py <action> [--path path] [--config] [configure_file]")
    actions = ["list_free_nodes", "wating_estimation", "print", "create", "submit", "array", "check_queue", "cancel_job", 
            "job_info", "print_avalable_steps"]
    parser = argparse.ArgumentParser()
    parser.add_argument("action", nargs="+", choices=actions, help="One or more actions to perform.")
    parser.add_argument("-p", "--path", type=str, help="Path to the working directory or resource.")
    parser.add_argument("-c", "--config", type=str, help="Path to a configuration file.")
    parser.add_argument("--id", type=str, help="Job ID (for details or cancellation).")
    parser.add_argument("-s", "--steps", help="Use steps (1,2,3 or 1-3")
    parser.add_argument("-a", "--array", action="store_true", help="Run job array based on step and path")
    parser.add_argument("-d", "--dependency_step", help="A step that blocks rest of steps")
#    check_action(parser.parse_args().action)
    return parser.parse_args()


#args = parser.parse_args()
#part_values = parse_part_values(args.part)
