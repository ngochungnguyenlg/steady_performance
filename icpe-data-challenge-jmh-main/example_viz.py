import re
import sys
import random
import json
import pandas as pd
import matplotlib.pyplot as plt
from glob import glob
from os.path import basename
import os
import numpy as np
from ultis import *
import ruptures as rpt
from crops import pelt_crops, segmentations


figs = "./figs"
class Benchmark:
    """Measurements from a single benchmark"""

    def __init__(self, filename):
        self.filename = filename
        # print(filename)

        m = re.match(r'(?P<org>[^_]+)__(?P<proj>[^#]+)#'
                     r'(?P<method>[^#]+)#(?P<params>.*)\.json',
                     basename(filename))
        if m is None:
            print('Unable to parse: {}'.format(filename), file=sys.stderr)
            return

        for k, v in m.groupdict().items():
            # print("k = {}, v = {}".format(k, v))
            setattr(self, k, v)

    def get_repository(self):
        return '{}/{}'.format(self.org, self.proj)

    def get_repository_url(self):
        return 'https://github.com/{}'.format(self.get_repository())

    def set_revision(self, revisions_csv):
        df = pd.read_csv(revisions_csv)
        df = df[df['repository'] == self.get_repository()]
        self.revision = df['tag'].iloc[0]

    def get_measurements(self):
        with open(self.filename) as f:
            self.measurements = json.load(f)
        return self.measurements
    
    def point_out_steady(self, fork_idx, visualize = False, min_p=4, max_p=1000, min_size=100, jump=100,max_iterator=200, method='l2'):
        fork_data= np.array(self.get_measurements()[fork_idx])*1000
        fil_fork_data = filter_outliers(fork_data, window_size=200)
        fil_fork_data=fil_fork_data.flatten()        
        segmentations = crops_algorithm(fil_fork_data, min_p, max_p, method=method, min_size=min_size, jump=jump, max_iterator=max_iterator)
        change_points = []
        n_change_points = []
        costs = []
        algo = rpt.Pelt(model=method).fit(fil_fork_data)
        total_point = 30
        penalties = [*segmentations]
        penalties.remove(max_p)
        
        if len(segmentations) <= total_point:
            penalties = np.linspace(min_p, max(penalties), total_point-len(segmentations))
        
        for penalty in penalties:
            result = algo.predict(pen=penalty)
            change_points += result
            n_change_points.append(len(result)-1)
            costs.append(calculate_segs_cost(result, penalty=penalty, model=algo))
            
        combined = list(zip(penalties, n_change_points, costs))
        combined_sorted = sorted(combined, key=lambda x: x[0])
                
        sorted_penalties, sorted_n_change_points, costs = zip(*combined_sorted)
        
        sorted_penalties = np.array(sorted_penalties)
        costs = np.array(costs)
        sorted_n_change_points = np.array(sorted_n_change_points)
        elbow = find_elbow(sorted_penalties, sorted_n_change_points) #penalty
        if elbow == None:
            elbow = penalties[0]
                    
        numper_change_point_optimize = len(algo.predict(pen=elbow))-1
        #using elbow to find best segments
        algo = rpt.Pelt(model=method, min_size=200).fit(fil_fork_data)
        segs = algo.predict(pen=elbow)
        stable = find_stable_segment(fil_fork_data, segs)
        
        if visualize:
            if not os.path.exists(figs):
                os.makedirs(figs)
            plt.figure(figsize=(10, 6))
            if elbow!=None:
                plt.scatter(numper_change_point_optimize,elbow, marker="*", color="red")
                plt.axvline(x=numper_change_point_optimize, color='r', linestyle='--')
            plt.plot(sorted_n_change_points, penalties, marker='o')
            plt.xlabel('Number of Changepoints')
            plt.ylabel('Penalty')
            plt.title('Penalty vs Number of Changepoints')
            plt.savefig(figs + "/"+ str(fork_idx)+"_fig5.png", dpi=300)                  
            plt.close()
            plt.plot(fil_fork_data)
            
            same_points = []
            for beta, segments in segmentations.items():
                for cp in segments:
                    if cp not in same_points:
                        same_points.append(cp)
                        plt.axvline(x=cp, color='r', linestyle='--')
            
            plt.savefig(figs + "/"+ str(fork_idx)+"_fig1.png", dpi=300)
            plt.close()
        return stable
    
    def steady_percentage(self, num):
        def cal_per_fork(stable_list, length):
            stable_length = 0
            for seg in stable_list:
                stable_length += seg[1]-seg[0]+1
            start_point = stable_list[-1][0]
            print(stable_length/length,  start_point)
            return stable_length/length, start_point
        data  = self.get_measurements()
        
        store_stable_percent= []
        store_start_stable_point= []
        for i in range(num):
            length = len(data[i])
            stable_range = self.point_out_steady(i, visualize=True)
            stable_percentage, start_stable_point = cal_per_fork(stable_range, length)
            store_stable_percent.append(stable_percentage)
            store_start_stable_point.append(start_stable_point)
        
        average_percent =  np.mean(store_stable_percent) 
        average_start_point=  np.mean(store_start_stable_point) 
        return average_percent, average_start_point
        

        
def get_benchmarks(data_dir):
    """Get all the jsons we can find and try to parse their names"""
    return [Benchmark(f) for f in glob('{}/*.json'.format(data_dir))]

def plot(benchmark):
    """Simple plot just to check that everything is in place"""

    # Get the measurements
    data = benchmark.get_measurements()

    # Create as many axes as the number of forks in the data
    fig, axs = plt.subplots(len(data), figsize=(15, 3 * len(data)))

    # Put some info about the benchmark in the title
    title = 'project {} at rev {}\nmethod: {}'.format(
            benchmark.get_repository(), benchmark.revision, benchmark.method)
    if benchmark.params != '':
        title += '\nparams: {}'.format(benchmark.params)
    fig.suptitle(title)

    # Make a scatter plot for each fork
    for i, fork in enumerate(data):
        x = [x for x in range(0, len(fork))]
        axs[i].scatter(x, fork)
        axs[i].set_title('Fork: {}'.format(i))
        axs[i].set_xlabel('iteration')
        axs[i].set_ylabel('Average exec time (seconds)')

    # Save
    plt.tight_layout(rect=[0, 0, 1, 0.98])
    plotfile = basename(benchmark.filename).replace('.json', '.png')
    plt.savefig(plotfile, dpi=300)
    print('Plot saved to: {}'.format(plotfile))

def random_viz(data_dir, revisions):
    """Randomly pick a benchmark to visualize"""

    # Find the measurements data
    benchmarks = get_benchmarks(DATA_DIR)

    # Randomly select a benchmark to visualize
    bench = random.choice(benchmarks)

    # Associate a revision to the benchmark
    bench.set_revision(REVISIONS)

    # Plot it
    plot(bench)
    
def fig_5(directory: str, file_name: tuple, randome_draw: bool, fork_n=2, min_size=2, jum=5 ):
    org = file_name[0]
    proj = file_name[1]
    method = params = None
    if len(file_name)==3:      
        method = file_name[2]
    elif len(file_name)==4:
        params = file_name[3]
        method = file_name[2]
    
    matched_files = []
        
    pattern = re.compile(rf'{org}__{proj}.*\.json')

    # Duyệt qua tất cả các file trong thư mục
    for filename in os.listdir(directory):
        # Kiểm tra xem file có khớp với mẫu không
        if pattern.match(filename) and method == None and params == None:
            matched_files.append(filename)
        elif pattern.match(filename) and method != None and params == None:
            if filename.find(method) != -1:
                matched_files.append(filename)
        elif pattern.match(filename) and method != None and params != None:
            if filename.find(method) != -1 and filename.find(params) != -1:
                matched_files.append(filename)
        else:
            pass
    
    def draw(filename,):
        ben = Benchmark(filename=directory+"/"+filename)
        # ben.point_out_steady(4, visualize=True, min_size=20)
        result = ben.steady_percentage(10)
        print(result)
        exit(0)
        data = ben.get_measurements()
        fork_data = np.array(data[fork_n])
        print(fork_data)
        model = "l2"  # Model cost function
        penalty = 0  # Giá trị penalty khởi tạo (có thể điều chỉnh)
        
        algo = rpt.Pelt(model=model).fit(fork_data)
        penalties = np.linspace(5, 160, 30)
        num_changepoints = []

        performances = []        
        for penalty in penalties:
            algo = rpt.Pelt(model="rbf", min_size=min_size, jump=jum).fit(fork_data)
            result = algo.predict(pen=penalty)
            print(len(result)-1)
            num_changepoints.append(len(result) - 1)
            performances += result

        plt.figure(figsize=(10, 6))
        plt.plot(num_changepoints, penalties, marker='o')
        plt.xlabel('Number of Changepoints')
        plt.ylabel('Penalty')
        plt.title('Penalty vs Number of Changepoints')
        plt.grid(True)
        plt.show()
        plt.close()
        
        plt.plot(fork_data)
        for cp in performances:
            plt.axvline(x=cp, color='r', linestyle='--')
        plt.xlabel('Iterator')
        plt.ylabel('Performance (Time)')
        plt.title('Performance Over Iterations')
        plt.grid(True)
        plt.legend()
        plt.show()
        
    if len(matched_files) and randome_draw:
        selected_file = random.choices(matched_files)
        draw(selected_file)
    elif len(matched_files) and not randome_draw:
        for fn in matched_files:
            draw(fn)
    return matched_files
    
if __name__ == "__main__":

    DATA_DIR = 'icpe-data-challenge-jmh-main/timeseries'
    REVISIONS = 'icpe-data-challenge-jmh-main/benchmarks_revision.csv'

    # random_viz(DATA_DIR, REVISIONS)
    #eclipse__eclipse-collections#org.eclipse.collections.impl.jmh.map.ChainMapPutTest.ec
    fig_5(directory=DATA_DIR, file_name = ("eclipse", "eclipse-collections", "org.eclipse.collections.impl.jmh.map.ChainMapPutTest.ec", "isPresized=true"), randome_draw=False)
