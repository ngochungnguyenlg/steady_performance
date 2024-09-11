import numpy as np
from icpe.config import *
import ruptures as rpt

def estimate_warmup_time(fork_m, w, wi):
    #this algorithm is from the paper to detect warmup time
    #algorithm 1.
    wt = 0
    sit = 0
    sic = 0
    for e in fork_m:
        ni = int(e/100)
        it = ni*e
        sit += it
        if sit >w:
            wt += sit
            sit = 0
            sic += 1
        if sic >= wi:
            break
    return wt

def select_performance_measurements(wt, r, i, M):
    #This algorithm is from the paper to detect M^{conf}.
    #algorithm 2.
    sit = 0
    sic = 0
    t = 0
    Mconf = []
    for e in M:
        ni = e / 100 
        it = ni * e
        if t >= wt:
            sit += it
            Mconf.append(e)
            if sit > r:
                sit = 0
                sic += 1
            if sic >= i:
                break
        t += e
    return Mconf


def filter_outliers(data, window_size=200):
    data = np.array(data)
    window_size = min(window_size, len(data))  
    noops = np.ceil(0.1 / np.mean(data))
    skip = int(np.ceil(window_size / np.log1p(noops)))

    indexes = list(range(skip, len(data)))
    outliers = set()

    def sliding_window(seq, window_size):
        for i in range(len(seq) - window_size + 1):
            yield seq[i:i + window_size]

    def create_isoutlier(series, window):
        series_ = series[window]
        median = np.median(series_)
        p1, p99 = [np.quantile(series_, q) for q in [0.01, 0.99]]
        threshold = 3 * (p99 - p1)

        def _isoutlier(i):
            value = series[i]
            return abs(value - median) > threshold

        return _isoutlier

    for window in sliding_window(indexes, window_size):
        isoutlier = create_isoutlier(data, window)
        outliers.update([i for i in window if isoutlier(i)])

    filtered_data = np.setdiff1d(np.arange(len(data)), list(outliers))
    
    return data[filtered_data]

def calculate_segs_cost(tau, penalty, model):
    tau = [0] + tau
    total_cost = 0
    
    for i in range(1, len(tau)):
        start = tau[i-1]+1
        end = tau[i]        
        if start < end-1:
            try:
                cost = model.cost.error(start, end-1) + penalty
                total_cost += cost 
            except Exception as e:
                print("Error happen, because of {}, start at {}, end at {}".format(e, start, end-1))
    return total_cost

def crops_algorithm(data, beta_min, beta_max, method = "l2", jump=5, min_size=2, max_iterator = 1000):
    '''
    params:
    
    '''
    if len(data) == 0:
        return {}
    model = rpt.Pelt(model=method, min_size=min_size, jump=jump).fit(data)
    def run_pelt(penalty):
        return model.predict(pen=penalty)    
            
    beta_set = {(beta_min, beta_max)}
    try:
        segmentations = {beta_min: run_pelt(beta_min), beta_max: run_pelt(beta_max)}
    except Exception as e:
        print (e, data)
    n= max_iterator
    while beta_set and n>0:
        beta0, beta1 = beta_set.pop()
        try:
            segs_0 = run_pelt(beta0)
            segs_1 = run_pelt(beta1)
            m0 = len(segs_0)
            m1 = len(segs_1)
            if m0 >m1 + 1:
                beta_int = (calculate_segs_cost(segs_1,beta1, model) - calculate_segs_cost(segs_0,beta0, model)) /(m0-m1)
                
                segments_int = run_pelt(beta_int)
                
                if len(segments_int) != m1:
                    beta_set.update([(beta0, beta_int), (beta_int, beta1)])
                    segmentations[beta_int] = segments_int
                print(len(beta_set), beta_int, m0, m1, beta0, beta1)
        except Exception as e:
            print (e, data)
        n -=1
    return segmentations


def find_elbow(penalties, costs):
    from kneed import KneeLocator
    kneedle = KneeLocator(penalties, costs, curve='convex', direction='decreasing')
    elbow_point = kneedle.elbow
    return elbow_point


def calculate_mean_and_ci(data):
    """Calculate the mean and 95% confidence interval of the mean."""
    mean = np.mean(data)
    ci = 1.96 * np.std(data, ddof=1) / np.sqrt(len(data))
    return mean, ci

def is_within_tolerance(last_mean, last_ci, new_mean, new_ci, tolerance=0.05):
    """Check if mean2 is within the tolerance of mean1."""
    last_low = last_mean-last_ci
    last_up = last_mean+last_ci
    
    new_low = new_mean-new_ci
    new_up = new_mean+new_ci
    
    low_difer = abs(last_low-new_low)/last_low
    up_difer = abs(last_up-new_up)/last_up
    
    return low_difer<tolerance and up_difer < tolerance

def find_stable_segment(data, segments):
    """Find stable segments based on the 5% tolerance on the CI of the mean."""
    stable_segments = []
    unstable_segments = []
    segments = sorted(segments)

    print("--->",segments, len(data))
    if len(segments) <2:
        return [(0, segments[0])], segments[0] >500
    elif len(segments)==2:
        return [(segments[0], segments[1])], abs(segments[1]-segments[0])>500
    if config["using_last_stable_segmet"]==False:
        max_length_seg = -float("inf")
        max_length_seg_idx = -1
        #ignore the first segnment
        for idx in range(0,len(segments)-1):
            if abs(segments[idx+1] - segments[idx]) > max_length_seg:
                max_length_seg = abs(segments[idx+1] - segments[idx])
                max_length_seg_idx = idx
        mean, last_ci = calculate_mean_and_ci(data[segments[max_length_seg_idx]:segments[max_length_seg_idx+1]])
        stable_segments.append((segments[max_length_seg_idx], segments[max_length_seg_idx+1]))   
    else:
        mean, last_ci = calculate_mean_and_ci(data[segments[-2]:segments[-1]])

    
    end_stabel_range = 0
    for i in range(len(segments) - 2, -1, -1):
        segment_data = data[segments[i]:segments[i+1]]
        last_mean, ci = calculate_mean_and_ci(segment_data)
        
        if is_within_tolerance(mean, ci, last_mean,last_ci):
            if (segments[i], segments[i+1]) not in stable_segments:
                stable_segments.append((segments[i], segments[i+1]))
                if end_stabel_range < segments[i+1]:
                    end_stabel_range = segments[i+1]
        else:
            unstable_segments.append((segments[i], segments[i+1]))
    is_steady = True
    print("--->", abs(segments[-1] - segments[-2]), abs(segments[-1] - segments[-2]) <= 500 and config["using_last_stable_segmet"])
    if abs(segments[-1] - segments[-2]) <= 500 and config["using_last_stable_segmet"]:
        is_steady=False
    elif  not config["using_last_stable_segmet"] and abs(segments[max_length_seg_idx+1]-segments[max_length_seg_idx])<500:
        is_steady=False

    stable_segments = sorted(stable_segments, key=lambda x: x[0], reverse=True)
    return stable_segments, is_steady

if __name__ == "__main__":
    np.random.seed(0)
    N = np.random.randint(100,200,5)
    noise = np.random.uniform(-3,3,5)
    data = np.concatenate([np.random.normal(loc+noise[loc], 0.5, size=N[loc]) for loc in range(5)], axis=0)
    data = data.flatten()
    beta_min = 4
    beta_max = 10^5
    segmentations = crops_algorithm(data, beta_min, beta_max, max_iterator=10, method="l1")

    import matplotlib.pyplot as plt
    plt.plot(data)
    for beta, segments in segmentations.items():
        for cp in segments:
            plt.axvline(x=cp, color='r', linestyle='--')
    plt.show()