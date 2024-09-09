import multiprocessing.pool
import re
import sys
import random
import json
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')  # Use non-GUI backend
from glob import glob
from os.path import basename
import os
import numpy as np
from icpe.ultis import *
import ruptures as rpt
from multiprocessing import Process,  Queue, Manager, Lock
import multiprocessing
from typing import List, Tuple
from functools import partial
import math
from icpe.config import config
from icpe.crops import *

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

        self.name = ""
        for k, v in m.groupdict().items():
            # print("k = {}, v = {}".format(k, v))
            setattr(self, k, v)
            self.name += k+"_"+v
        self.lock =  Lock()

    def get_name(self):
        return self.name
    
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
    
    def point_out_steady(self, fork_idx, visualize = config["visualize"], min_p=4, max_p=1000, min_size=100, jump=100, max_iterator=200, method='l2'):
        mean = np.mean(self.get_measurements()[fork_idx])
        exponent_10 = math.floor(math.log10(mean))
        fork_data= np.array(self.get_measurements()[fork_idx])*(10/(10**exponent_10))
        fil_fork_data = filter_outliers(fork_data, window_size=200)
        fil_fork_data=fil_fork_data.flatten()
        if len(fil_fork_data) == 0:
            return [], 0, 0
        # segmentations = crops_algorithm(fil_fork_data, min_p, max_p, method=method, min_size=min_size, jump=jump, max_iterator=max_iterator)
        segmentations = pelt_crops(data=fil_fork_data, beta_min=min_p, beta_max=max_p, max_iterations=max_iterator, jump=jump, min_size=min_size, method=method)
        change_points = []
        n_change_points = []
        costs = []
        algo = rpt.Pelt(model=method, min_size=min_size, jump=jump).fit(fil_fork_data)
        total_point = 30
        
        if visualize == True:
            penalties = [*segmentations]
            if max_p in penalties:
                penalties.remove(max_p)
            
            penalties = list(np.unique(np.array(penalties)))
            if len(segmentations) <= total_point and len(penalties)>1:
                penalties = np.linspace(min_p, max(penalties), total_point-len(penalties))
            
            for penalty in penalties:
                result = algo.predict(pen=penalty)
                change_points += result
                n_change_points.append(len(result)-1)
                costs.append(calculate_segs_cost(result, penalty=penalty, model=algo))
        else:
            penalties = [*segmentations]
            n_change_points = []
            for pen in penalties:
                n_change_points.append(len(segmentations[pen]))
            
        combined = list(zip(penalties, n_change_points))
        combined_sorted = sorted(combined, key=lambda x: x[0])
                
        sorted_penalties, sorted_n_change_points = zip(*combined_sorted)
        
        sorted_penalties = np.array(sorted_penalties)
        sorted_n_change_points = np.array(sorted_n_change_points)
        elbow = find_elbow(sorted_penalties, sorted_n_change_points) #penalty
        if elbow == None:
            elbow = penalties[0]
                    
        numper_change_point_optimize = len(algo.predict(pen=elbow))-1
        #using elbow to find best segments
        segs = algo.predict(pen=elbow)
        stable = find_stable_segment(fil_fork_data, segs)
        
        total_seg = len(segs)-1
        stable_seg_n = len(stable)
        isconsitent = True
        if total_seg - stable_seg_n >1:
            isconsitent = False
        if visualize:
            self.lock.acquire()
            if not os.path.exists(figs):
                os.makedirs(figs)
            plt.figure(figsize=(10, 6))
            if elbow!=None:
                try:
                    plt.scatter(numper_change_point_optimize,elbow, marker="*", color="red")
                    plt.axvline(x=numper_change_point_optimize, color='r', linestyle='--')
                except:
                    print("---------->", numper_change_point_optimize, elbow, type(numper_change_point_optimize), type(elbow))
            plt.plot(sorted_n_change_points, penalties, marker='o')
            plt.xlabel('Number of Changepoints')
            plt.ylabel('Penalty')
            plt.title('Penalty vs Number of Changepoints')
            plt.savefig(figs + "/"+ str(fork_idx)+"_fig5.png", dpi=300)                  
            plt.close()
            plt.plot(fil_fork_data)
            
            for cp in segs:
                plt.axvline(x=cp, color='r', linestyle='--')
            
            plt.savefig(figs + "/"+ str(fork_idx)+"_fig1.png", dpi=300)
            plt.close()
            self.lock.release()
        return stable, isconsitent, stable_seg_n > 1
    
    def steady_percentage(self, num):
        def cal_per_fork(stable_list, length):
            stable_length = 0
            for seg in stable_list:
                stable_length += seg[1] - seg[0] + 1
            start_point = stable_list[-1][0]
            print(stable_length/length,  start_point)
            return stable_length/length, start_point
        data  = self.get_measurements()
        store_stable_percent= []
        store_start_stable_point= []
        store_is_consistent = []
        store_is_steady = []
        for i in range(num):
            length = len(data[i])
            stable_range, is_consistent, steady = self.point_out_steady(i, visualize=True, min_size=200, max_p=10**5, jump=5, max_iterator=20)
            stable_percentage, start_stable_point = cal_per_fork(stable_range, length)
            store_stable_percent.append(stable_percentage)
            store_start_stable_point.append(start_stable_point)
            store_is_consistent.append(is_consistent)
            store_is_steady.append(steady)
        
        average_percent =  np.mean(store_stable_percent) 
        average_start_point=  np.mean(store_start_stable_point) 
        return average_percent, average_start_point, is_consistent, steady

    def process_task(self, index: int, fidx) -> Tuple[float, int]:
        print("analysing file {} thread {} starting".format(fidx,index))
        data = self.get_measurements()
        length = len(data[index])
        stable_range, is_consistent, steady = self.point_out_steady(index, visualize=config["visualize"], min_size=200, max_p=10**5, jump=5, max_iterator=20)
        stable_percentage, start_stable_point = self.cal_per_fork(stable_range, length)
        print("analysing file {} thread {} end".format(fidx,index))
        return stable_percentage, sum(data[index][0:start_stable_point]), is_consistent, steady

    def cal_per_fork(self, stable_list, length):
        stable_length = 0
        for seg in stable_list:
            stable_length += seg[1] - seg[0] + 1
        start_point = stable_list[-1][0]
        print(stable_length / length, start_point)
        return stable_length / length, start_point 
        
    def steady_percentage_by_mul_thread(self, num, file_idx):
        store_stable_percent = []
        store_start_stable_point = []
        store_is_consistent = []
        store_is_steady = []

        process_task_with_file_idx = partial(self.process_task, fidx=file_idx)
        with  multiprocessing.pool.ThreadPool(processes=num) as pool:
            results = pool.map(process_task_with_file_idx, range(num))
        for stable_percentage, start_stable_point, is_consistent, is_steady in results:
            store_stable_percent.append(stable_percentage)
            store_start_stable_point.append(start_stable_point)
            store_is_consistent.append(is_consistent)
            store_is_steady.append(is_steady)

        # average_percent = np.mean(store_stable_percent)
        # average_start_point = np.mean(store_start_stable_point)
        return store_stable_percent, store_start_stable_point, store_is_consistent, store_is_steady

        
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
    benchmarks = get_benchmarks(data_dir)

    # Randomly select a benchmark to visualize
    bench = random.choice(benchmarks)

    # Associate a revision to the benchmark
    bench.set_revision(revisions)

    # Plot it
    plot(bench)
    
def analyze_benchmark(ben, idx):
    print("analysis at idx {} of {}".format(idx, ben.get_name()))
    return ben.get_name(), ben.steady_percentage(10)

def analyze_benchmark_thread(ben, idx, queue):
    print("analysis at idx {} of {}".format(idx, ben.get_name()))
    queue.put([ben.get_name(), ben.steady_percentage_by_mul_thread(10, idx)])

def exe_process(fr: int, to: int, benchmarks: List, queue):
    processes = []
    for idx in range(fr, to):
        p = Process(target=analyze_benchmark_thread, args=(benchmarks[idx], idx, queue))
        processes.append(p)
        p.start()
    for p in processes:
        p.join()

def main_processes(benchmarks: List, max_worker: int):
    result = {}
    manager = Manager()
    queue = manager.Queue()

    for i in range(0, len(benchmarks), max_worker):
        exe_process(i, min(i + max_worker, len(benchmarks)), benchmarks, queue)
        while not queue.empty():
            name, data = queue.get()
            if isinstance(data, Exception):
                print(f'Benchmark {name} generated an exception: {data}')
            else:
                result[name] = data
                print(f"Analysis complete for {name}")
        try:
            df = pd.DataFrame(result.items(), columns=['Benchmark', 'Steady Percentage'])
            df.to_csv(config['output_dir'], index=False)
        except Exception as e:
            print("Error during save csv file {}".format(e))
            
def analysis_data(data_dir, resvision, action, num_worker, range_data=[]):
    print(
        data_dir,
        resvision,
        action,
        num_worker,
        range_data
    )
    if action == "example":
        random_viz(data_dir, resvision)
    elif action == "debug":
        benchmarks = get_benchmarks(data_dir)[range_data[0]:range_data[1]]
        benchmarks[0].steady_percentage(1)
    else:
        if range_data!=None:
            benchmarks = get_benchmarks(data_dir)[range_data[0]:range_data[1]]
        else:
            benchmarks = get_benchmarks(data_dir)
        main_processes(benchmarks, num_worker)
    
            
if __name__ == "__main__":
    DATA_DIR = './icpe/timeseries'
    REVISIONS = './icpe/benchmarks_revision.csv'
    analysis_data(DATA_DIR,REVISIONS, "run", 50)