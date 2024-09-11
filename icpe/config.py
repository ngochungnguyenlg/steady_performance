
config = {
    "using_last_stable_segmet": True,
    "output_dir" : "./csv/output.csv",
    "figs" : "figs",
    "visualize": True
}
import os
if not os.path.exists(config['figs']):
    os.makedirs(config['figs'])
if not os.path.exists(config['output_dir'].split("/")[1]):
    os.makedirs(config['output_dir'].split("/")[1])