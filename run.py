import argparse
from icpe.analysis_csv import *
from icpe.example_viz import *

def parser_input():
    parser = argparse.ArgumentParser(description="The input for experiment of icpe.")
    parser.add_argument(
        '-i', '--input', 
        type=str, 
        required=True, 
        help='Running option: analysis data (data), analysis result (result)'
    )
    
    parser.add_argument(
        '--dlink', 
        type=str, 
        default="./icpe/timeseries", 
    )
    
    parser.add_argument(
        '--rlink', 
        type=str, 
        nargs='+',
        default="./icpe/benchmarks_revision.csv",
    )
    
    parser.add_argument(
        '--range', 
        type=int, 
        nargs='+'
    )
    
    parser.add_argument(
        '--nworker', 
        type=int, 
        default= 10,
        help='Number of worker for analysis'
    )
    
    parser.add_argument(
        '--action', 
        type=str, 
        default= "run",
    )

    parser.add_argument(
        '-o', '--option', 
        choices=['fig5', 'fig6', 'fig8'], 
        help='Choose an option from predefined choices'
    )
    return parser.parse_args()  

def run_data_analysis(args):
    analysis_data(
        data_dir=args.dlink,
        resvision=args.rlink[0],
        action = args.action,
        num_worker=args.nworker,
        range_data=args.range
        )

def run_result_analysis(args):
    
    if args.action == "rq1":
        research_question_1(args.rlink[0], args.dlink)
    elif args.action == "rq2":
        research_question_2(args.rlink[0], args.dlink)
    elif args.action == "rq3":
        research_question_3(args.rlink[0], args.dlink)
    elif args.action == "rq4":
        research_question_4(args.rlink[0], args.dlink)
    elif args.action == "rq5":
        research_question_5(args.rlink[0], args.dlink)
    elif args.action == "rq6":
        research_question_6(args.rlink, args.dlink)
    else:
        print("wrong selection !!!")

def main():
    args = parser_input()
    if args.input == "data":
        run_data_analysis(args)
    elif args.input == "result":
        run_result_analysis(args)
    else:
        print("wrong selection !!!\n select 'data' for data analysis,\n 'result' for result analysis")
if __name__ == "__main__":
    main()
