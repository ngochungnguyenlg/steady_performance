import pandas as pd
import os
import ast
import matplotlib.pyplot as plt
# plt.style.use('ggplot')
import seaborn as sns
import pandas as pd
import numpy as np
import scienceplots
from icpe.config import config
# matplotlib.rcParams['font.size']=12.5
plt.style.use(["science", "ieee"])
import matplotlib
matplotlib.rcParams.update({'font.size': 12})

class survey_result:
    def __init__(self, csv_link = "result.csv", raw_data_link = "timeseries"):
        self.df = pd.read_csv(csv_link)
        self.methods = self.get_all_method(raw_data_link)
        self.methods.sort(key=str.lower)
        
    def get_all_method(self, link):
        file_name = os.listdir(link)
        methods = [m[m.find("_")+2:m.find("#")].replace(" ","") for m in file_name]
        methods = list(set(methods))
        return methods
    
    def detach_rq1_result_cell(self, cell):
        cell = ast.literal_eval(cell)
        #steady percent, start steady point, is_consistent, is_stable (has aleast a steady seg)
        return list(cell[0]), list(cell[1]), list(cell[2]), list(cell[3])
    
    def draw_table_heat_map(self, data, name):
        df = pd.DataFrame(data, index=self.methods)
        plt.figure(figsize=(6, 10))  
        ax = sns.heatmap(df, 
                        annot=True, 
                        cmap='Greys', 
                        fmt='.2f', 
                        linewidths=0.5, 
                        linecolor='black', 
                        cbar_kws={"shrink": .8}, 
                        vmin=df.min().min(),      
                        vmax=df.max().max()
                        ) 
        plt.xticks(rotation=0, ha='center')           
        plt.tight_layout()  
        plt.savefig(name, dpi=300)
        # plt.show()
        plt.close()

    def draw_f8(self, data):
        fig = plt.figure(figsize =(10, 7))
        ax = fig.add_subplot(111)
        bp = ax.boxplot(data)
        ax.set_xticklabels(self.methods, rotation=90)
        plt.yscale('log')
        plt.ylim(0.1, None)
        plt.ylabel("Time to reach the steady state (sec), (log scale)")
        plt.savefig("./figs/fig.8.png", dpi=300)
        plt.close()
        
    def question_number_one(self):
        df = self.df.copy()
        results = {}
        box_plot_start_point =[]
        
        heat_map_table_data ={"inconsistent": [], "steady state":[]}
        heat_map_table_data_6a ={"no steady state": [], "steady state":[]}
        for method in self.methods:
            method_result = df[df["Benchmark"].str.contains(method)]
            measurement_result = []
            benchmark_inconsistent = 0 # when some fork cannot reach steady state, result no 4
            start_point = []
            no_steady_fork = 0
            for cel in method_result["Steady Percentage"]:
                steady_percentage, steady_point, is_consistent, has_steady =  self.detach_rq1_result_cell(cel)
                if False in is_consistent:
                    benchmark_inconsistent += 1
                for r in is_consistent:
                    if r==False:
                        no_steady_fork += 1
                start_point+= steady_point
            if len(method_result["Steady Percentage"]) >0:
                heat_map_table_data["inconsistent"].append(100.0*benchmark_inconsistent/len(method_result["Steady Percentage"]))
                heat_map_table_data["steady state"].append(100.0-100.0*benchmark_inconsistent/len(method_result["Steady Percentage"]))
                heat_map_table_data_6a["no steady state"].append(100*no_steady_fork/(10*len(method_result["Steady Percentage"])))
                heat_map_table_data_6a["steady state"].append(100-100*no_steady_fork/(10*len(method_result["Steady Percentage"])))
            else:
                heat_map_table_data["inconsistent"].append(0.0)
                heat_map_table_data["steady state"].append(0.0)
                heat_map_table_data_6a["no steady state"].append(0.0)
                heat_map_table_data_6a["steady state"].append(0.0)
                
            results[method] = measurement_result
            box_plot_start_point.append(np.array(start_point))
        
        self.draw_table_heat_map(heat_map_table_data, "./figs/fig6b.png")
        self.draw_table_heat_map(heat_map_table_data_6a,"./figs/fig6a.png")
        self.draw_f8(box_plot_start_point)

def research_question_1(raw_data_link, csv_link=config["output_dir"]):
    s = survey_result(csv_link, raw_data_link = raw_data_link)
    s.question_number_one()
        
def research_question_2():
    raise NotImplementedError("This method is not implemented")

def research_question_3():
    raise NotImplementedError("This method is not implemented")

def research_question_4():
    raise NotImplementedError("This method is not implemented")

def research_question_5():
    raise NotImplementedError("This method is not implemented")     
        
if __name__ == "__main__":
    s = survey_result(csv_link="./result.csv", raw_data_link="./icpe-data-challenge-jmh-main/timeseries")
    s.question_number_one()
            